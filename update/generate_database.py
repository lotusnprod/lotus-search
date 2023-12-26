#!/usr/bin/env python3
import logging
import pickle
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    database = {}

    with open(path / "database_chemo.pkl", "rb") as f:
        database.update(pickle.load(f))
    logging.info("Finished integrating chemo")

    with open(path / "database_taxo.pkl", "rb") as f:
        database.update(pickle.load(f))
    logging.info("Finished integrating taxo")

    with open(path / "database_biblio.pkl", "rb") as f:
        database.update(pickle.load(f))
    logging.info("Finished integrating biblio")

    with open(path / "database.pkl", "wb") as f:
        pickle.dump(database, f)
    logging.info("Finished dumping")


if __name__ == "__main__":
    run(Path("data"))
