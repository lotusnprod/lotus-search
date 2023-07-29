#!/usr/bin/env python3

from pathlib import Path
from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path) -> None:
    query = """
    SELECT DISTINCT ?structure ?structure_smiles ?canonical_smiles WHERE {
        ?structure wdt:P703 ?taxon;
                   wdt:P233 ?canonical_smiles.
        OPTIONAL {
            ?structure wdt:P2017 ?structure_smiles.
        }
    }
    """

    t = wd_sparql_to_csv(query)

    with open(root / "smiles.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
