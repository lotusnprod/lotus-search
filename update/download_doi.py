#!/usr/bin/env python3

from pathlib import Path

from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def run(root: Path, retry: int = 3) -> None:
    query = """
    SELECT DISTINCT ?reference ?reference_doi WITH {
      SELECT DISTINCT ?reference WHERE {
        # We use InChIKey (P235) instead of SMILES as some of them are incomplete.
        ?structure wdt:P235 [].
        ?structure p:P703/(prov:wasDerivedFrom/pr:P248) ?reference.
      }
      # LIMIT 1000 # Test
    } AS %refs WHERE {
      SELECT DISTINCT ?reference ?reference_doi WHERE { 
        INCLUDE %refs
                ?reference wdt:P356 ?reference_doi.
      }
    }
    """

    t = wd_sparql_to_csv(query)
    if "java.util.concurrent.TimeoutException" in t:
        if retry > 0:
            print("  Failed to download DOIs, retrying...")
            run(root, retry - 1)
            return
        else:
            raise TimeoutError("Failed to download DOIs tum tum tum....")

    with open(root / "dois.csv", "w") as f:
        f.write(remove_wd_entity_prefix(t))
