#!/usr/bin/env python3

from pathlib import Path

from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path, retry: int = 3) -> None:
    query = """
    SELECT DISTINCT ?structure ?structure_smiles ?canonical_smiles WHERE {
        # Using InChIKey (P235) to recognize chemicals. 
        # Could also be
        # P31 wd:Q113145171 `type of a chemical entity`
        # P31 wd:Q59199015 `group of stereoisomers`
        ?structure wdt:P235 [].
        ?structure wdt:P703 ?taxon.
          
        # All P2017 should also have P233 but some of them are not complete.
        OPTIONAL {
            ?structure wdt:P233 ?canonical_smiles.
        }           
        OPTIONAL {
            ?structure wdt:P2017 ?structure_smiles.
        }
    }
    """

    t = wd_sparql_to_csv(query)
    if "java.util.concurrent.TimeoutException" in t:
        if retry > 0:
            print("  Failed to download smiles, retrying...")
            run(root, retry - 1)
            return
        else:
            raise TimeoutError("Failed to download smiles tum tum tum....")

    with open(root / "smiles.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
