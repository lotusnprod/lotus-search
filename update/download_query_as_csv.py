#!/usr/bin/env python3
import logging
from pathlib import Path

from update.common import WD_URL, remove_wd_entity_prefix, sparql_to_csv
import sys


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def run(path: Path, query_file: Path, output_file: Path,
        retry: int = 2, url: str = WD_URL) -> None:
    with open(query_file, "r") as qf:
        query = qf.read()

    t = sparql_to_csv(query=query, url=url)
    if "java.util.concurrent.TimeoutException" in t:
        if retry > -1:
            logging.warning(f"Failed to download query {query_file}, retrying...")
            run(path, query_file, output_file, retry - 1, url)
            return
        else:
            if url != WD_URL:
                logging.warning(f"Timeout for query from file {query_file} on url {url} retrying with {WD_URL}")
                run(path, query_file, output_file, 2, WD_URL)
                return
            logging.error(f"Timeout for query from file {query_file}")
            raise TimeoutError(f"Failed to download query {query_file}....")

    with open(path / output_file, "w") as f:
        f.write(remove_wd_entity_prefix(text=t))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise ValueError(f"Usage: python3 {sys.argv[0]} <query_file> <output_file>")
    run(path=Path("data"), query_file=Path(sys.argv[1]), output_file=Path(sys.argv[2]))