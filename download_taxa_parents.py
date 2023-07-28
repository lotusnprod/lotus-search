#!/usr/bin/env python3

import requests
import time

start = time.time()

url = "https://query.wikidata.org/sparql"
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

r = requests.get(url, params={'query': query},
                 headers={'Accept': 'text/csv'})

with open("data/taxa.csv", "w") as f:
    f.write(r.text.replace("http://www.wikidata.org/entity/Q", ""))

print(f"Time elapsed: {time.time() - start:.2f}s")
