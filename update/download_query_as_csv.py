import logging
import sys
from pathlib import Path
from time import sleep

from update.common import (
    QLEVER_URL,
    WD_URL,
    remove_wd_entity_prefix_and_Q,
    sparql_to_csv,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def run(
    path: Path,
    query_file: Path,
    output_file: Path,
    retry: int = 3,
    url: str = WD_URL,
) -> None:
    with open(query_file) as qf:
        query = qf.read()

    t = sparql_to_csv(
        query=query,
        url=url,
        as_post=(retry == 1),
    )  # The last retry we do a POST just in case the cache was polluted
    if "java.util.concurrent.TimeoutException" in t:
        if retry > 1:
            logging.warning(f"Failed to download query {query_file}, retrying...")
            sleep(5)
            run(path, query_file, output_file, retry - 1, url)
            return
        else:
            logging.warning(
                f"Timeout for query from file {query_file} on url {url} retrying with {QLEVER_URL}",
            )
            run(path, query_file, output_file, 1, QLEVER_URL)
            return
        logging.error(f"Timeout for query from file {query_file}")
        raise TimeoutError(f"Failed to download query {query_file}....")

    logging.info("Keeping unique lines")
    lines = t.splitlines()
    header = lines[0]  # Save the header line
    unique_content_lines = list(set(lines[1:]))  # Get unique content lines

    # Process each line separately
    processed_lines = [header]  # Start with the header
    for line in unique_content_lines:
        processed_line = remove_wd_entity_prefix_and_Q(text=line)
        if processed_line is not None:  # Only add non-None results
            processed_lines.append(processed_line)

    # Join all processed lines
    t_processed = "\n".join(processed_lines)

    logging.info(f"Writing to {path / output_file}")
    with open(path / output_file, "w") as f:
        f.write(t_processed)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise ValueError(f"Usage: python3 {sys.argv[0]} <query_file> <output_file>")
    run(path=Path("data"), query_file=Path(sys.argv[1]), output_file=Path(sys.argv[2]))
