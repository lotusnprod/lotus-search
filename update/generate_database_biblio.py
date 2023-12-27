#!/usr/bin/env python3
import logging
import pickle
from csv import DictReader
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    reference_doi: dict[int, str] = {}

    with open(path / "references.csv", "r") as f:
        reader = DictReader(f)
        dict_doi: dict = {int(row["reference"]): row["reference_doi"] for row in reader}

    database = {
        "reference_doi": dict_doi,
    }

    with open(path / "database_biblio.pkl", "wb") as f:
        pickle.dump(database, f)


if __name__ == "__main__":
    run(Path("data"))
