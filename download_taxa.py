#!/usr/bin/env python3

import requests
import time

start = time.time()

url = "https://query.wikidata.org/sparql"
query = """
SELECT DISTINCT ?taxon ?taxon_name WHERE {
  ?compound wdt:P703 ?taxon.         # Found in taxon
  ?taxon wdt:P225 ?taxon_name.       # Get scientific name of the taxon
}
"""

r = requests.get(url, params={'query': query},
                 headers={'Accept': 'text/csv'})

with open("data/taxa.csv", "w") as f:
    f.write(r.text.replace("http://www.wikidata.org/entity/Q", ""))

print(f"Time elapsed: {time.time() - start:.2f}s")
