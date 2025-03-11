from requests import request

WD_URL = "https://query.wikidata.org/sparql"
QLEVER_URL = "https://qlever.cs.uni-freiburg.de/api/wikidata"


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
    """
    Removes the Wikidata entity prefix and then, if the identifier starts with a 'Q',
    removes the 'Q' so that only the numeric portion remains.
    
    If the identifier starts with an 'L' (a lexeme), returns None to signal that
    lexemes should be excluded.
    """
    cleaned = remove_wd_entity_prefix(text)
    if cleaned.startswith("L"):
        # Exclude lexeme identifiers.
        return None
    if cleaned.startswith("Q"):
        # Remove the leading "Q" for Q-entities.
        cleaned = cleaned[1:]
    return cleaned
