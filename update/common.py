from requests import request

# Allows to use queries without SERVICE <>
WD_URL = "https://query-legacy-full.wikidata.org/sparql"
QLEVER_URL = "https://qlever.dev/api/wikidata"


def sparql_to_csv(query: str, url: str = WD_URL, as_post: bool = False) -> str:
    method = "POST" if as_post else "GET"
    return request(
        method,
        url,
        params={"query": query},
        headers={
            "Accept": "text/csv",
            "Accept-Encoding": "gzip,deflate",
            "User-Agent": "LOTUS project database dumper",
        },
    ).text


def remove_wd_entity_prefix(text: str) -> str:
    return text.replace("http://www.wikidata.org/entity/", "")


def remove_wd_entity_prefix_and_Q(text: str) -> str:
    cleaned = remove_wd_entity_prefix(text)
    if cleaned.startswith("L"):  # Exclude lexemes
        return None
    return cleaned.replace("Q", "")
