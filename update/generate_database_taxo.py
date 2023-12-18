#!/usr/bin/env python3
import logging
import pickle
from csv import DictReader
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


# See https://www.wikidata.org/wiki/Q2576881
def convert_to_int_safe(s: str):
    try:
        result = int(s)
        return result
    except ValueError:
        logging.warning(f"{s} is not a valid integer.")
        return None


def run(path: Path) -> None:
    taxon_direct_parents: dict[int, set[int]] = {}
    taxon_names: dict[int, str] = {}
    taxon_ranks: dict[int, set[int]] = {}
    taxon_children: dict[int, set[int]] = {}
    taxon_parents_with_distance: dict[int, dict[int, int]] = {}
    ranks_names: dict[int, str] = {}

    with open(path / "taxa.csv", "r") as t:
        reader = DictReader(t)
        for row in reader:
            taxon_name = row["taxon_name"]
            taxon_id = convert_to_int_safe(row["taxon"])
            parent_id = convert_to_int_safe(row["parent"])
            taxon_rank_id = convert_to_int_safe(row["taxon_rank"])

            if parent_id not in taxon_children:
                taxon_children[parent_id] = set()
            if taxon_id not in taxon_direct_parents:
                taxon_direct_parents[taxon_id] = set()
            if taxon_id not in taxon_names:
                taxon_names[taxon_id] = taxon_name
            if taxon_id not in taxon_ranks:
                taxon_ranks[taxon_id] = set()
            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = {}

            taxon_direct_parents[taxon_id].add(parent_id)
            taxon_children[parent_id].add(taxon_id)
            taxon_parents_with_distance[taxon_id][parent_id] = 1
            taxon_ranks[taxon_id].add(taxon_rank_id)

        logging.info(f" Found {len(taxon_direct_parents)} taxa")

    with open(path / "taxa_parents.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            taxon_id = convert_to_int_safe(row["taxon"])
            taxon_name = str(row["taxon_name"])
            taxon_rank_id = convert_to_int_safe(row["taxon_rank"])
            relative_id = convert_to_int_safe(row["relative"])
            relative_name = str(row["relative_name"])
            relative_rank_id = convert_to_int_safe(row["relative_rank"])
            distance = convert_to_int_safe(row["distance"])

            if relative_id not in taxon_children:
                taxon_children[relative_id] = set()
            if taxon_id not in taxon_direct_parents:
                taxon_direct_parents[taxon_id] = set()
            if taxon_id not in taxon_names:
                taxon_names[taxon_id] = taxon_name
            if relative_id not in taxon_names:
                taxon_names[relative_id] = relative_name
            if taxon_id not in taxon_ranks:
                taxon_ranks[taxon_id] = set()
            if relative_id not in taxon_ranks:
                taxon_ranks[relative_id] = set()
            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = {}

            if distance == 1:
                taxon_direct_parents[taxon_id].add(relative_id)

            taxon_children[relative_id].add(taxon_id)
            # We also add the children of the ones above
            if taxon_id in taxon_children:
                for child_id in taxon_children[taxon_id]:
                    taxon_children[relative_id].add(child_id)

            taxon_ranks[taxon_id].add(taxon_rank_id)

            taxon_ranks[relative_id].add(relative_rank_id)

            taxon_parents_with_distance[taxon_id][relative_id] = distance

    with open(path / "taxa_ranks.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            ranks_names[int(row["rank"])] = row["rankLabel"]

    with open(path / "taxa_all.csv", "r") as f:
        reader = DictReader(f)
        dict_all_taxa: dict = {row["taxon"]: row["taxon_name"] for row in reader}

    database = {
        "taxonomy_direct_parents": taxon_direct_parents,
        "taxonomy_names": dict_all_taxa,
        "taxonomy_ranks": taxon_ranks,
        "taxonomy_children": taxon_children,
        "taxonomy_parents_with_distance": taxon_parents_with_distance,
        "taxonomy_ranks_names": ranks_names,
    }

    with open(path / "database_taxo.pkl", "wb") as f:
        pickle.dump(database, f)


if __name__ == "__main__":
    run(Path("data"))
