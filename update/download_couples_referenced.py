#!/usr/bin/env python3

from pathlib import Path

from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path, retry: int = 3) -> None:
    query = """
    SELECT DISTINCT ?compound ?taxon ?reference WHERE {
      ?compound p:P703 ?statement;
        # We use InChIKey (P235) instead of SMILES as some of them are incomplete.
        wdt:P235 [].
      ?statement ps:P703 ?taxon;
        (prov:wasDerivedFrom/pr:P248) ?reference.
    }
    """

    t = wd_sparql_to_csv(query)
    if "java.util.concurrent.TimeoutException" in t:
        if retry > 0:
            print("  Failed to download taxonomy step 1, retrying...")
            run(root, retry - 1)
            return
        else:
            raise TimeoutError("Failed to download taxonomy step 1 tum tum tum....")

    with open(root / "couples.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
