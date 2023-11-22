#!/usr/bin/env python3
import csv
import pickle
from pathlib import Path


def run(root: Path) -> None:
    database = {}

    with open(root / "database_chemo.pkl", "rb") as f:
        database.update(pickle.load(f))
    print("Finished integrating chemo")

    with open(root / "database_taxo.pkl", "rb") as f:
        database.update(pickle.load(f))
    print("Finished integrating taxo")

    # TODO references

    with open(root / "database.pkl", "wb") as f:
        pickle.dump(database, f)
    print("Finished dumping")


if __name__ == "__main__":
    run(Path("data"))
