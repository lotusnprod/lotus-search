from pathlib import Path
from requests import get

WD_URL = "https://query.wikidata.org/sparql"
QLEVER_URL = "https://qlever.cs.uni-freiburg.de/api/wikidata"


def wd_sparql_to_csv(query: str, url: str = WD_URL) -> str:
    return get(
        url, params={"query": query}, headers={"Accept": "text/csv"}
    ).text


def remove_wd_entity_prefix(text: str) -> str:
    return text.replace("http://www.wikidata.org/entity/Q", "")


def run_query_to_csv(query_file: str, output_file: str, root: Path = Path("data"), retry: int = 3, url: str = WD_URL) -> None:
    with open(query_file, "r") as qf:
        query = qf.read()

    t = wd_sparql_to_csv(query=query, url=url)
    if "java.util.concurrent.TimeoutException" in t:
        if retry > 0:
            print("  Failed to download query, retrying...")
            run_query_to_csv(query_file, output_file, root, retry - 1, url)
            return
        else:
            raise TimeoutError("Failed to download query tum tum tum....")

    with open(root / output_file, "w") as f:
        f.write(remove_wd_entity_prefix(text=t))
