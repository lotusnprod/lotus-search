#!/usr/bin/env python3

from pathlib import Path
from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path) -> None:
    query = """
    SELECT DISTINCT ?taxon ?taxon_name ?parent_taxon ?parent_taxon_name WITH {
      SELECT DISTINCT ?taxon WHERE {
        ?compound wdt:P703 ?taxon.
      }
    } AS %taxa WITH {
      SELECT ?taxon ?taxon_name ?parent_taxon WHERE {
        INCLUDE %taxa
        ?taxon wdt:P225 ?taxon_name;
               wdt:P171* ?parent_taxon.
      }
    } AS %taxa_2 WHERE {
      SELECT ?taxon ?taxon_name ?parent_taxon ?parent_taxon_name WHERE {
        INCLUDE %taxa_2
        ?parent_taxon wdt:P225 ?parent_taxon_name;
      }
    }
    """

    t = wd_sparql_to_csv(query)

    with open(root / "taxa.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
