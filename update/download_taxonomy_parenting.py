import csv
import pickle
from io import StringIO
from pathlib import Path

from update.common import remove_wd_entity_prefix, wd_sparql_to_csv

query_taxa = """
SELECT DISTINCT ?taxon ?taxon_name ?taxon_rank ?parent WHERE {
    ?compound wdt:P703 ?taxon.
    ?taxon wdt:P171 ?parent;
           wdt:P105 ?taxon_rank;
           wdt:P225 ?taxon_name.
  }
"""

query_final = """
SELECT ?taxon ?taxon_name ?taxon_rank ?relative ?relative_name ?relative_rank ?distance WITH {
  SELECT DISTINCT ?taxon WHERE {
    ?compound wdt:P703 ?subtax.
    ?subtax wdt:P171 ?taxon.
  }
} AS %taxa WHERE {
  SELECT DISTINCT ?taxon ?taxon_name ?taxon_rank ?relative ?relative_name ?relative_rank (COUNT(DISTINCT ?rel) AS ?distance) WHERE  { 
  INCLUDE %taxa
  ?taxon wdt:P171* ?rel;
         wdt:P105 ?taxon_rank;
         wdt:P225 ?taxon_name.
  ?rel wdt:P171+ ?relative.
  ?relative wdt:P105 ?relative_rank;
            wdt:P225 ?relative_name.
  } GROUP BY ?taxon ?taxon_name ?taxon_rank ?relative ?relative_name ?relative_rank
}
"""

query_ranks = """
SELECT ?rank ?rankLabel WHERE {
    ?rank wdt:P31 wd:Q427626.
    SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
"""


def run(root: Path) -> None:
    t = remove_wd_entity_prefix(wd_sparql_to_csv(query_taxa))

    list_of_couples = [x.strip().split(",") for x in t.split("\n")[1:] if x != ""]
    taxon_direct_parents = {}
    taxon_names = {}
    taxon_ranks = {}
    taxon_children = {}
    taxon_parents_with_distance = {}
    ranks_names = {}

    for couple in list_of_couples:
        taxon_id, taxon_name, taxon_rank_id, parent_id = couple
        taxon_id = int(taxon_id)
        parent_id = int(parent_id)
        taxon_rank_id = int(taxon_rank_id)

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

    print(f" Found {len(taxon_direct_parents)} taxa")

    t = remove_wd_entity_prefix(wd_sparql_to_csv(query_final))
    parents = {}

    reader = csv.reader(StringIO(t))
    reader.__next__()
    for line in reader:
        taxon_id, taxon_name, taxon_rank_id, relative_id, relative_name, relative_rank, distance = line
        taxon_id = int(taxon_id)
        relative_id = int(relative_id)
        taxon_rank_id = int(taxon_rank_id)
        distance = int(distance)

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

        taxon_ranks[relative_id].add(relative_rank)

        taxon_parents_with_distance[taxon_id][relative_id] = distance

    t = remove_wd_entity_prefix(wd_sparql_to_csv(query_ranks))
    reader = csv.reader(StringIO(t))
    reader.__next__()
    for line in reader:
        ranks_names[int(line[0])] = line[1]

    database = {
        "taxonomy_direct_parents": taxon_direct_parents,
        "taxonomy_names": taxon_names,
        "taxonomy_ranks": taxon_ranks,
        "taxonomy_children": taxon_children,
        "taxonomy_parents_with_distance": taxon_parents_with_distance,
        "taxonomy_ranks_names": ranks_names
    }

    with open(root / "database_taxo.pkl", "wb") as f:
        pickle.dump(database, f)
