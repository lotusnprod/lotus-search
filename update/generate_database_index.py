#!/usr/bin/env python3
import csv
import logging
from pathlib import Path

from storage.storage import Storage
from update.taxo_helper import convert_to_int_safe, generate_taxon_parents_with_distance

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    storage = Storage(path)
    storage.drop_and_create_tables()
    taxa = set()  # So we can get only the parenting that we care about
    triplets = []
    dedup = set()
    with open(path / "triplets.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        ref_index = headers.index("reference")
        struct_index = headers.index("structure")
        taxon_index = headers.index("taxon")

        for row in reader:
            ref = int(row[ref_index])
            struct = int(row[struct_index])
            taxon = int(row[taxon_index])
            taxa.add(taxon)
            if (ref, struct, taxon) in dedup:
                continue
            dedup.add((ref, struct, taxon))
            triplets.append(
                {
                    "reference_id": int(row[ref_index]),
                    "structure_id": int(row[struct_index]),
                    "taxon_id": int(row[taxon_index]),
                }
            )
    logging.info(" Processed triplets")

    references_dict = {}
    journals_dict = {}
    with open(path / "references.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        ref_index = headers.index("reference")
        doi_index = headers.index("reference_doi")
        title_index = headers.index("reference_title")
        date_index = headers.index("reference_date")
        journal_index = headers.index("reference_journal")
        journal_title_index = headers.index("journal_title")

        for row in reader:
            references_dict[int(row[ref_index])] = {
                "doi": row[doi_index],
                "title": row[title_index],
                "date": row[date_index],
                "journal": row[journal_index],
            }
            journals_dict[int(row[journal_index])] = row[journal_title_index]

    logging.info("Processed references and journals")

    with open(path / "structures_table.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        # Get the indices of the columns
        id_index = headers.index("structure")
        smiles_index = headers.index("structure_smiles")
        # Add indices for new columns
        smiles_no_stereo_index = headers.index("structure_smiles_no_stereo")
        inchi_index = headers.index("structure_inchi")
        inchi_no_stereo_index = headers.index("structure_inchi_no_stereo")
        inchikey_index = headers.index("structure_inchikey")
        inchikey_no_stereo_index = headers.index("structure_inchikey_no_stereo")
        formula_index = headers.index("structure_formula")

        structures_dict = {}
        for row in reader:
            struct_id = int(row[id_index])
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

    logging.info(" Processed structures")
    with open(path / "taxa_names.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        taxon_index = headers.index("taxon")
        name_index = headers.index("taxon_name")

        taxo_names_dict = {int(row[taxon_index]): row[name_index] for row in reader}

    logging.info(" Processed taxa names")

    # Eventually TODO add taxa_names_com

    taxon_ranks_dict = {}

    with open(path / "ranks_names.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        rank_index = headers.index("rank")
        label_index = headers.index("rankLabel")

        ranks_names = [
            {"id": int(row[rank_index]), "name": row[label_index]} for row in reader
        ]

    logging.info(" Processed rank names")

    with open(path / "taxa_ranks.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        taxon_index = headers.index("taxon")
        rank_index = headers.index("taxon_rank")

        for row in reader:
            rank_value = convert_to_int_safe(row[rank_index])
            if rank_value is not None:
                taxon_ranks_dict[int(row[taxon_index])] = {rank_value}

    logging.info(" Processed taxa ranks")
    taxo_ranks = []
    for taxon, ranks in taxon_ranks_dict.items():
        for rank in ranks:
            taxo_ranks.append({"id": taxon, "rank_id": rank})
    taxo_names = []
    for taxon, name in taxo_names_dict.items():
        taxo_names.append({"id": taxon, "name": name})

    structures = list(structures_dict.values())

    references = []
    for ref, values in references_dict.items():
        references.append(
            {
                "id": ref,
                "doi": values["doi"],
                "title": values["title"],
                "date": values["date"],
                "journal": values["journal"],
            }
        )

    journals = []
    for journal, title in journals_dict.items():
        journals.append({"id": journal, "title": title})

    logging.info(" Processed dicts")
    storage.upsert_taxo_parenting(generate_taxon_parents_with_distance(path))
    logging.info(" Taxo parenting inserted")

    storage.upsert_triplets(triplets)
    logging.info(" Triplets inserted")
    storage.upsert_structures(structures)
    logging.info(" Structures inserted")
    storage.upsert_references(references)
    logging.info(" References inserted")
    storage.upsert_journals(journals)
    logging.info(" Journals inserted")
    storage.upsert_taxo_names(taxo_names)
    logging.info(" Taxo names inserted")
    storage.upsert_rank_names(ranks_names)
    logging.info(" Rank names inserted")
    storage.upsert_taxo_ranks(taxo_ranks)
    logging.info(" Taxo ranks inserted")
    logging.info("Finished generating index database")


if __name__ == "__main__":
    run(Path("data"))
