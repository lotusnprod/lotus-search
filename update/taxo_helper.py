#!/usr/bin/env python3
import csv
import logging
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


def generate_taxon_parents_with_distance(path: Path) -> dict[int, dict[int, int]]:
    taxon_parents_with_distance: dict[int, dict[int, int]] = {}

    with open(path / "taxa.csv", "r") as t:
        reader = csv.reader(t)
        headers = next(reader)
        taxon_index = headers.index("taxon")
        parent_index = headers.index("parent")
        for row in reader:
            taxon_id = convert_to_int_safe(row[taxon_index])
            parent_id = convert_to_int_safe(row[parent_index])

            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = dict()

            taxon_parents_with_distance[taxon_id][parent_id] = 1

    with open(path / "taxa_parents.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        taxon_index = headers.index("taxon")
        relative_index = headers.index("relative")
        distance_index = headers.index("distance")

        for row in reader:
            taxon_id = convert_to_int_safe(row[taxon_index])
            relative_id = convert_to_int_safe(row[relative_index])
            distance = convert_to_int_safe(row[distance_index])
            if taxon_id is None or relative_id is None or distance is None:
                continue

            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = dict()

            taxon_parents_with_distance[taxon_id][relative_id] = distance
        logging.info(f" Found {len(taxon_parents_with_distance)} taxa with parents.")

    return taxon_parents_with_distance
