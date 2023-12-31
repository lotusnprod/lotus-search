#!/usr/bin/env python3
import logging
import pickle
from collections import defaultdict
from csv import DictReader
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


# See https://www.wikidata.org/wiki/Q2576881
def convert_to_int_safe(s: str) -> int | None:
    try:
        result = int(s)
        return result
    except ValueError:
        logging.error(f"{s} is not a valid integer.")
        return None


def run(path: Path) -> None:
    taxon_ranks: dict[int, set[int]] = defaultdict(set)
    taxon_parents_with_distance: dict[int, dict[int, int]] = {}

    with open(path / "taxa.csv", "r") as t:
        reader = DictReader(t)
        for row in reader:
            taxon_id = convert_to_int_safe(row["taxon"])
            parent_id = convert_to_int_safe(row["parent"])
            taxon_rank_id = convert_to_int_safe(row["taxon_rank"])

            if taxon_id not in taxon_ranks:
                taxon_ranks[taxon_id] = set()
            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = dict()

            taxon_parents_with_distance[taxon_id][parent_id] = 1
            taxon_ranks[taxon_id].add(taxon_rank_id)


    with open(path / "taxa_parents.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            taxon_id = convert_to_int_safe(row["taxon"])
            relative_id = convert_to_int_safe(row["relative"])
            distance = convert_to_int_safe(row["distance"])

            if taxon_id is None or relative_id is None or distance is None:
                continue

            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = dict()

            taxon_parents_with_distance[taxon_id][relative_id] = distance
        logging.info(f" Found {len(taxon_parents_with_distance)} taxa with parents.")

    database = {
        "taxonomy_parents_with_distance": taxon_parents_with_distance,
    }

    # TODO make a csv instead of a pkl and probably just do all that in the index task
    with open(path / "database_taxo.pkl", "wb") as f:
        pickle.dump(database, f)


if __name__ == "__main__":
    run(Path("data"))
