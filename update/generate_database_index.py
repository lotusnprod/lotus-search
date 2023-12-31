#!/usr/bin/env python3
import logging
import pickle
from csv import DictReader
from pathlib import Path

from storage.storage import Storage
from update.generate_database_taxo import convert_to_int_safe

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    storage = Storage(path)
    triplets = []
    with open(path / "triplets.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            r = int(row["reference"])
            s = int(row["structure"])
            t = int(row["taxon"])
            triplets.append({"reference_id": r, "structure_id": s, "taxon_id": t})

    storage.upsert_triplets(triplets)

    with open(path / "references.csv", "r") as f:
        reader = DictReader(f)
        references = [
            {"id": int(row["reference"]), "doi": row["reference_doi"]} for row in reader
        ]

    with open(path / "smiles_processed.csv", "r") as f:
        reader = DictReader(f)
        structures = [
            {"id": int(row["structure"]), "smiles": row["structure_smiles"]}
            for row in reader
        ]

    with open(path / "taxa_names.csv", "r") as f:
        reader = DictReader(f)
        taxo_names = [
            {"id": int(row["taxon"]), "name": row["taxon_name"]} for row in reader
        ]
    dict_taxon_ranks = {}
    with open(path / "taxa.csv", "r") as t:
        reader = DictReader(t)
        for row in reader:
            taxon_id = int(row["taxon"])
            taxo_names += [{"id": taxon_id, "name": row["taxon_name"]}]
            dict_taxon_ranks[taxon_id] = {int(row["taxon_rank"])}

    with open(path / "ranks_names.csv", "r") as f:
        reader = DictReader(f)
        ranks_names = [
            {"id": int(row["rank"]), "name": row["rankLabel"]} for row in reader
        ]

    with open(path / "database_taxo.pkl", "rb") as f:
        database_taxo = pickle.load(f)
        storage.upsert_taxo_parenting(database_taxo["taxonomy_parents_with_distance"])

    with open(path / "taxa_ranks.csv", "r") as f:
        reader = DictReader(f)

        for row in reader:
            rank_value = convert_to_int_safe(row["taxon_rank"])
            if rank_value is not None:
                dict_taxon_ranks[int(row["taxon"])] = set([rank_value])
        logging.info(f" Found {len(dict_taxon_ranks)} taxa with rank.")

    taxo_ranks = []
    for taxon, ranks in dict_taxon_ranks.items():
        for rank in ranks:
            taxo_ranks.append({"id": taxon, "rank_id": rank})

    storage.upsert_structures(structures)
    logging.info(" Structures inserted")
    storage.upsert_references(references)
    logging.info(" References inserted")
    storage.upsert_taxo_names(taxo_names)
    logging.info(" Taxo names inserted")
    storage.upsert_rank_names(ranks_names)
    logging.info(" Rank names inserted")
    storage.upsert_taxo_ranks(taxo_ranks)
    logging.info(" Taxo ranks inserted")
    logging.info("Finished generating index database")


if __name__ == "__main__":
    run(Path("data"))
