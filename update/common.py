import requests

WD_URL = "https://query.wikidata.org/sparql"


def wd_sparql_to_csv(query: str) -> str:
    return requests.get(
        WD_URL, params={"query": query}, headers={"Accept": "text/csv"}
    ).text


def remove_wd_entity_prefix(text: str) -> str:
    return text.replace("http://www.wikidata.org/entity/Q", "")
