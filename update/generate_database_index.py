#!/usr/bin/env python3
import logging
from csv import DictReader
from pathlib import Path

from storage.storage import Storage

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

    storage.upsert_structures(structures)
    storage.upsert_references(references)
    logging.info("Finished generating index database")


if __name__ == "__main__":
    run(Path("data"))
