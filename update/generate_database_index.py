import csv
import logging
from pathlib import Path

from storage.storage import Storage
from update.taxo_helper import generate_taxon_parents_with_distance

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def _read_taxa_names(path: Path) -> dict[int, str]:
    taxo_names_path = path / "taxa_names.csv"
    if not taxo_names_path.exists():
        logging.error("taxa_names.csv not found at %s", taxo_names_path)
        return {}
    good: dict[int, str] = {}
    skipped = 0
    with open(taxo_names_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        try:
            headers = next(reader)
        except StopIteration:
            logging.warning("taxa_names.csv empty")
            return {}
        try:
            taxon_index = headers.index("taxon")
            name_index = headers.index("taxon_name")
        except ValueError:
            logging.error("Missing required columns in taxa_names.csv: %s", headers)
            return {}
        for row in reader:
            if len(row) <= max(taxon_index, name_index):
                skipped += 1
                continue
            raw_id = row[taxon_index].strip()
            name = row[name_index]
            # Fast digit check then fallback to try/except
            if not raw_id or not raw_id.isdigit():
                try:
                    tax_id = int(raw_id)
                except Exception:
                    skipped += 1
                    continue
            else:
                tax_id = int(raw_id)
            good[tax_id] = name
    if skipped:
        logging.info(
            "taxa_names.csv: %d entries parsed, %d rows skipped (invalid id)",
            len(good),
            skipped,
        )
    return good


def _read_references_and_journals(path: Path) -> tuple[list[dict], list[dict]]:
    """Parse references.csv into references and journals lists safely.

    Skips malformed / short rows instead of raising IndexError / ValueError.
    Preserves last occurrence for duplicate reference IDs (prior behavior).
    """
    ref_path = path / "references.csv"
    if not ref_path.exists():
        logging.error("references.csv not found at %s", ref_path)
        return [], []

    references_dict: dict[int, dict] = {}
    journals_dict: dict[int, str] = {}
    skipped = 0

    with open(ref_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        try:
            headers = next(reader)
        except StopIteration:
            logging.warning("references.csv empty")
            return [], []

        required = [
            "reference",
            "reference_doi",
            "reference_title",
            "reference_date",
            "reference_journal",
            "journal_title",
        ]
        missing = [h for h in required if h not in headers]
        if missing:
            logging.error("Missing required columns in references.csv: %s", missing)
            return [], []

        idx = {name: headers.index(name) for name in required}
        max_required_index = max(idx.values())

        for row in reader:
            # Skip short or empty rows
            if len(row) <= max_required_index:
                skipped += 1
                continue
            try:
                ref_id_raw = row[idx["reference"]].strip()
                journal_id_raw = row[idx["reference_journal"]].strip()
                if not ref_id_raw.isdigit() or not journal_id_raw.isdigit():
                    skipped += 1
                    continue
                ref_id = int(ref_id_raw)
                journal_id = int(journal_id_raw)
                references_dict[ref_id] = {
                    "id": ref_id,
                    "doi": row[idx["reference_doi"]],
                    "title": row[idx["reference_title"]],
                    "date": row[idx["reference_date"]],
                    "journal": journal_id,
                }
                journals_dict[journal_id] = row[idx["journal_title"]]
            except Exception:  # broad: keep ingestion resilient
                skipped += 1
                continue

    if skipped:
        logging.info(
            "references.csv: %d references parsed, %d rows skipped (invalid / short)",
            len(references_dict),
            skipped,
        )

    references = list(references_dict.values())
    journals = [{"id": jid, "title": title} for jid, title in journals_dict.items()]
    return references, journals


def run(path: Path) -> None:
    storage = Storage(path)
    # TODO change this
    storage.drop_and_create_tables()
    taxa = set()  # So we can get only the parenting that we care about
    triplets = []
    dedup = set()
    with open(path / "triplets.csv") as f:
        reader = csv.reader(f)
        headers = next(reader)
        ref_index = headers.index("reference")
        struct_index = headers.index("structure")
        taxon_index = headers.index("taxon")

        for row in reader:
            # Skip malformed triplet row
            if len(row) <= max(ref_index, struct_index, taxon_index):
                continue
            try:
                ref = int(row[ref_index])
                struct = int(row[struct_index])
                taxon = int(row[taxon_index])
            except ValueError:
                continue
            taxa.add(taxon)
            if (ref, struct, taxon) in dedup:
                continue
            dedup.add((ref, struct, taxon))
            triplets.append({
                "reference_id": ref,
                "structure_id": struct,
                "taxon_id": taxon,
            })
    logging.info(" Processed triplets")

    # Rewritten robust parsing for references and journals
    references, journals = _read_references_and_journals(path)
    logging.info("Processed references and journals")

    structures_dict = {}
    with open(path / "structures_table.csv") as f:
        reader = csv.reader(f)
        headers = next(reader)
        id_index = headers.index("structure")
        smiles_index = headers.index("structure_smiles")
        smiles_no_stereo_index = headers.index("structure_smiles_no_stereo")
        inchi_index = headers.index("structure_inchi")
        inchi_no_stereo_index = headers.index("structure_inchi_no_stereo")
        inchikey_index = headers.index("structure_inchikey")
        inchikey_no_stereo_index = headers.index("structure_inchikey_no_stereo")
        formula_index = headers.index("structure_formula")
        max_struct_idx = max(
            id_index,
            smiles_index,
            smiles_no_stereo_index,
            inchi_index,
            inchi_no_stereo_index,
            inchikey_index,
            inchikey_no_stereo_index,
            formula_index,
        )
        for row in reader:
            if len(row) <= max_struct_idx:
                continue
            try:
                struct_id = int(row[id_index])
            except ValueError:
                continue
            structures_dict[struct_id] = {
                "id": struct_id,
                "smiles": row[smiles_index],
                "smiles_no_stereo": row[smiles_no_stereo_index],
                "inchi": row[inchi_index],
                "inchi_no_stereo": row[inchi_no_stereo_index],
                "inchikey": row[inchikey_index],
                "inchikey_no_stereo": row[inchikey_no_stereo_index],
                "formula": row[formula_index],
            }
    structures = list(structures_dict.values())
    logging.info(" Processed structures")

    descriptors_dict = {}
    with open(path / "descriptors_rdkit.csv") as f:
        reader = csv.reader(f)
        headers = next(reader)
        smiles_index = headers.index("smiles")
        descriptor_indices = range(1, len(headers))
        max_desc_idx = len(headers) - 1
        for row in reader:
            if len(row) <= max_desc_idx:
                continue
            smiles = row[smiles_index]
            if not smiles:
                continue
            struct_data = {}
            for i in descriptor_indices:
                value = row[i]
                struct_data[headers[i]] = float(value) if value else None
            struct_data["smiles"] = smiles
            descriptors_dict[smiles] = struct_data
    logging.info("Processed descriptors")

    taxo_names_dict = _read_taxa_names(path)
    taxo_names = [
        {"id": taxon, "name": name} for taxon, name in taxo_names_dict.items()
    ]
    logging.info(" Processed taxa names")

    ranks_names = []
    with open(path / "ranks_names.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                ranks_names.append({"id": int(row["rank"]), "name": row["rankLabel"]})
            except (ValueError, KeyError):
                continue
    logging.info(" Processed rank names")

    taxon_ranks_dict = {}
    with open(path / "taxa_ranks.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                taxon_id = int(row["taxon"])  # noqa: F841 (document structure for future multi-rank support)
                rank_value = int(row["taxon_rank"])  # noqa: F841
                taxon_ranks_dict[taxon_id] = {rank_value}
            except (ValueError, KeyError):
                continue
    taxo_ranks = [
        {"id": taxon, "rank_id": rank}
        for taxon, ranks in taxon_ranks_dict.items()
        for rank in ranks
    ]
    logging.info(" Processed taxa ranks")

    logging.info(" Generating taxonomy, might take long")
    storage.upsert_taxo_parenting(generate_taxon_parents_with_distance(path))
    logging.info(" Taxo parenting inserted")
    storage.upsert_triplets(triplets)
    logging.info(" Triplets inserted")
    storage.upsert_structures(structures)
    logging.info(" Structures inserted")
    logging.info(" Inserting descriptors, might take long")
    storage.upsert_structures_descriptors(descriptors_dict)
    logging.info(" Structures descriptors inserted")
    storage.upsert_references(references)
    logging.info(" References inserted")
    storage.upsert_journals(journals)
    logging.info(" Journals inserted")
    logging.info(" Inserting taxo names, might take long")
    storage.upsert_taxo_names(taxo_names)
    logging.info(" Taxo names inserted")
    logging.info(" Inserting taxo ranks, might take long")
    storage.upsert_rank_names(ranks_names)
    logging.info(" Rank names inserted")
    storage.upsert_taxo_ranks(taxo_ranks)
    logging.info(" Taxo ranks inserted")
    logging.info("Finished generating index database")


if __name__ == "__main__":
    run(Path("data"))
