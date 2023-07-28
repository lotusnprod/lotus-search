#!/usr/bin/env python3

import requests
import time

start = time.time()

url = "https://query.wikidata.org/sparql"
query = """
SELECT DISTINCT ?compound ?taxon ?reference WHERE {
  ?compound p:P703 ?statement;
    wdt:P233 ?canonical_smiles.
  ?statement ps:P703 ?taxon;
    (prov:wasDerivedFrom/pr:P248) ?reference.
}
"""

r = requests.get(url, params={'query': query},
                 headers={'Accept': 'text/csv'})

with open("data/couples.csv", "w") as f:
    f.write(r.text.replace("http://www.wikidata.org/entity/Q", ""))

print(f"Time elapsed: {time.time() - start:.2f}s")
