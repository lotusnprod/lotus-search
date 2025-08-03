#!/usr/bin/env python3
import csv
import logging
from pathlib import Path

from storage.storage import Storage
from update.taxo_helper import generate_taxon_parents_with_distance

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


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
            ref = int(row[ref_index])
            struct = int(row[struct_index])
            taxon = int(row[taxon_index])
            taxa.add(taxon)
            if (ref, struct, taxon) in dedup:
                continue
            dedup.add((ref, struct, taxon))
            triplets.append({
                "reference_id": int(row[ref_index]),
                "structure_id": int(row[struct_index]),
                "taxon_id": int(row[taxon_index]),
            })
    logging.info(" Processed triplets")

    references_dict = {}
    journals_dict = {}
    with open(path / "references.csv") as f:
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
    references = []
    for ref, values in references_dict.items():
        references.append({
            "id": ref,
            "doi": values["doi"],
            "title": values["title"],
            "date": values["date"],
            "journal": values["journal"],
        })

    journals = []
    for journal, title in journals_dict.items():
        journals.append({"id": journal, "title": title})
    logging.info("Processed references and journals")

    structures_dict = {}
    with open(path / "structures_table.csv") as f:
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
    structures = list(structures_dict.values())
    # structures = []
    # for struct, smiles in structures_dict.items():
    #     structures.append({"id": struct, "smiles": smiles})
    logging.info(" Processed structures")

    # TODO add all IDs and formatters (See #50)

    descriptors_dict = {}
    with open(path / "descriptors_rdkit.csv") as f:
        reader = csv.reader(f)
        headers = next(reader)
        smiles_index = headers.index("smiles")
        # Excluding the SMILES column
        descriptor_indices = range(1, len(headers))
        for row in reader:
            smiles = row[smiles_index]
            if smiles:
                # Create a dictionary for the current structure with SMILES as key
                struct_data = {}
                for i in descriptor_indices:
                    # Check if the value is not an empty string before converting to float
                    value = row[i]
                    if value:
                        struct_data[headers[i]] = float(value)
                    else:
                        # Handle empty value (perhaps set it to None or another default value)
                        logging.warning(f"Empty descriptor found in row: {row}")
                        struct_data[headers[i]] = None  # or any default value you prefer
                # Add SMILES as a separate key
                struct_data["smiles"] = smiles
                # Add the structure data to the descriptors dictionary
                descriptors_dict[smiles] = struct_data
            else:
                # Handle empty SMILES
                logging.warning(f"Empty SMILES found in row: {row}")
    logging.info("Processed descriptors")

    with open(path / "taxa_names.csv") as f:
        reader = csv.reader(f)
        headers = next(reader)
        taxon_index = headers.index("taxon")
        name_index = headers.index("taxon_name")
        taxo_names_dict = {int(row[taxon_index]): row[name_index] for row in reader}
    taxo_names = []
    for taxon, name in taxo_names_dict.items():
        taxo_names.append({"id": taxon, "name": name})
    logging.info(" Processed taxa names")

    # Eventually TODO add taxa_names_com

    # TODO add all IDs and formatters (See #50)

    # Process rank names
    ranks_names = []
    with open(path / "ranks_names.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                ranks_names.append({"id": int(row["rank"]), "name": row["rankLabel"]})
            except (ValueError, KeyError):
                logging.exception(f"Invalid row: {row}")
                continue

    # Process taxon names
    taxo_names = [{"id": taxon, "name": name} for taxon, name in taxo_names_dict.items()]
    logging.info(" Processed rank names")

    # Process taxon ranks
    taxon_ranks_dict = {}
    with open(path / "taxa_ranks.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                taxon_id = int(row["taxon"])
                rank_value = int(row["taxon_rank"])
                taxon_ranks_dict[taxon_id] = {rank_value}
            except (ValueError, KeyError):
                logging.exception(f"Invalid row: {row}")
                continue

    # Create final taxon ranks list
    taxo_ranks = [{"id": taxon, "rank_id": rank} for taxon, ranks in taxon_ranks_dict.items() for rank in ranks]

    logging.info(" Processed taxa ranks")
    logging.info(" Processed dicts")

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
