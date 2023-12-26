#!/usr/bin/env python3
import logging
import pickle
from csv import DictReader
from pathlib import Path

from storage.storage import Storage

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    storage = Storage(path)
    triplets = []
    with open(path / "couples.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            r = int(row["reference"])
            s = int(row["structure"])
            t = int(row["taxon"])
            triplets.append((r, s, t))

    storage.add_triplets(triplets)
    storage.close()
    logging.info("Finished generating index database")


if __name__ == "__main__":
    run(Path("data"))
