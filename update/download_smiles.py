#!/usr/bin/env python3

from pathlib import Path
from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path, retry: int = 3) -> None:
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
    if "java.util.concurrent.TimeoutException" in t:
        if retry > 0:
            print("  Failed to download smiles, retrying...")
            run(root, retry-1)
            return
        else:
            raise TimeoutError("Failed to download smiles tum tum tum....")


    with open(root / "smiles.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
