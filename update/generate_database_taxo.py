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
        logging.error(f"{s} is not a valid integer.")
        return None


def run(path: Path) -> None:
    taxon_direct_parents: dict[int, set[int]] = {}
    taxon_names: dict[int, str] = {}
    taxon_ranks: dict[int, set[int]] = {}
    taxon_children: dict[int, set[int]] = {}
    taxon_parents_with_distance: dict[int, dict[int, int]] = {}

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
        logging.info(f" Found {len(taxon_direct_parents)} taxa in which chemicals were found in.")

    with open(path / "taxa_names.csv", "r") as f:
        reader = DictReader(f)
        dict_taxon_names: dict = {int(row["taxon"]): row["taxon_name"] for row in reader}
        logging.info(f" Found {len(dict_taxon_names)} taxa with name.")

    with open(path / "taxa_ranks.csv", "r") as f:
        reader = DictReader(f)
        dict_taxon_ranks: dict = {int(row["taxon"]): convert_to_int_safe(row["taxon_rank"]) for row in reader}
        logging.info(f" Found {len(dict_taxon_ranks)} taxa with rank.")

    with open(path / "taxa_parents.csv", "r") as f:
        reader = DictReader(f)
        for row in reader:
            taxon_id = convert_to_int_safe(row["taxon"])
            relative_id = convert_to_int_safe(row["relative"])
            distance = convert_to_int_safe(row["distance"])

            if relative_id not in taxon_children:
                taxon_children[relative_id] = set()
            if taxon_id not in taxon_direct_parents:
                taxon_direct_parents[taxon_id] = set()
            if taxon_id not in dict_taxon_ranks:
                dict_taxon_ranks[taxon_id] = set()
            if relative_id not in dict_taxon_ranks:
                dict_taxon_ranks[relative_id] = set()
            if taxon_id not in taxon_parents_with_distance:
                taxon_parents_with_distance[taxon_id] = {}

            if distance == 1:
                taxon_direct_parents[taxon_id].add(relative_id)

            taxon_children[relative_id].add(taxon_id)
            # We also add the children of the ones above
            if taxon_id in taxon_children:
                for child_id in taxon_children[taxon_id]:
                    taxon_children[relative_id].add(child_id)

            taxon_parents_with_distance[taxon_id][relative_id] = distance
        logging.info(f" Found {len(taxon_parents_with_distance)} taxa with parents.")

    with open(path / "ranks_names.csv", "r") as f:
        reader = DictReader(f)
        ranks_names: dict = {int(row["rank"]): row["rankLabel"] for row in reader}
        logging.info(f" Found {len(ranks_names)} ranks with name.")

    database = {
        "taxonomy_direct_parents": taxon_direct_parents,
        "taxonomy_names": dict_taxon_names,
        "taxonomy_ranks": dict_taxon_ranks,
        "taxonomy_children": taxon_children,
        "taxonomy_parents_with_distance": taxon_parents_with_distance,
        "taxonomy_ranks_names": ranks_names,
    }

    with open(path / "database_taxo.pkl", "wb") as f:
        pickle.dump(database, f)


if __name__ == "__main__":
    run(Path("data"))
