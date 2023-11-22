import requests

WD_URL = "https://query.wikidata.org/sparql"
QLEVER_URL = "https://qlever.cs.uni-freiburg.de/api/wikidata"


def wd_sparql_to_csv(query: str, url: str = WD_URL) -> str:
    return requests.get(
        url, params={"query": query}, headers={"Accept": "text/csv"}
    ).text


def remove_wd_entity_prefix(text: str) -> str:
    return text.replace("http://www.wikidata.org/entity/Q", "")
