from requests import get, post

WD_URL = "https://query.wikidata.org/sparql"
QLEVER_URL = "https://qlever.cs.uni-freiburg.de/api/wikidata"


def sparql_to_csv(query: str, url: str = WD_URL, as_post: bool = False) -> str:
    if as_post:
        f = post
    else:
        f = get
    return f(
        url,
        params={"query": query},
        headers={
            "Accept": "text/csv",
            "Accept-Encoding": "gzip,deflate",
            "User-Agent": "LOTUS project database dumper",
        },
    ).text


def remove_wd_entity_prefix(text: str) -> str:
    return text.replace("http://www.wikidata.org/entity/Q", "")
