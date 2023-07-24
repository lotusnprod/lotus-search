#!/usr/bin/env python3

import requests
import time

start = time.time()

url = "https://query.wikidata.org/sparql"
query = """
SELECT DISTINCT ?structure ?structure_smiles ?canonical_smiles WHERE {
    ?structure wdt:P703 ?taxon;
               wdt:P233 ?canonical_smiles.
    OPTIONAL {
        ?structure wdt:P2017 ?structure_smiles.
    }
}
"""

r = requests.get(url, params={'query': query},
                 headers={'Accept': 'text/csv'})
with open("data/smiles.csv", "w") as f:
    f.write(r.text.replace("http://www.wikidata.org/entity/Q", ""))

print(f"Time elapsed: {time.time() - start:.2f}s")
