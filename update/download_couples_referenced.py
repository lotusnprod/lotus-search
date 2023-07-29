#!/usr/bin/env python3

from pathlib import Path
from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path) -> None:
    query = """
    SELECT DISTINCT ?compound ?taxon ?reference WHERE {
      ?compound p:P703 ?statement;
        wdt:P233 ?canonical_smiles.
      ?statement ps:P703 ?taxon;
        (prov:wasDerivedFrom/pr:P248) ?reference.
    }
    """

    t = wd_sparql_to_csv(query)

    with open(root / "couples.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
