import logging
import os
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

# Environment tunables (do not break tests; defaults are conservative)
MAX_RETRIES = int(os.getenv("LOTUS_UPDATE_MAX_RETRIES", "6"))
SWITCH_TO_QLEVER_AFTER = int(
    os.getenv("LOTUS_UPDATE_SWITCH_AFTER", "3"),
)  # attempts on WD before switching to POST and then QLever
SLEEP_TIME = 5  # seconds between retries (while still on original endpoint)
RATE_LIMIT_TOKENS = ["rate limit", "ratelimit", "too many requests"]


def _looks_rate_limited(payload: str) -> bool:
    if not payload:  # empty
        return True
    low = payload.lower()
    return any(tok in low for tok in RATE_LIMIT_TOKENS)


def _looks_timeout(payload: str) -> bool:
    return "java.util.concurrent.TimeoutException" in payload


def run(
    path: Path,
    query_file: Path,
    output_file: Path,
    retry: int = 3,  # retained for backwards compatibility (unused in new loop)
    url: str = WD_URL,
) -> None:
    """Download SPARQL results to CSV with order-preserving de-duplication.

        Retry logic:
          - Uses GET for the first SWITCH_TO_QLEVER_AFTER-1 attempts, then POST.
          - Switches from the original endpoint to QLever after SWITCH_TO_QLEVER_AFTER attempts.
          - Sleeps SLEEP_TIME seconds between retries while still on the original endpoint
            (i.e., only for attempt < SWITCH_TO_QLEVER_AFTER).

    Parameters
    ----------
    path : Path
        Path.
    query_file : Path
        Query file.
    output_file : Path
        Output file.
    retry : int
        Default is 3.
    url : str
        WD_URL. Default is WD_URL.
    """
    with open(query_file) as qf:
        query = qf.read()

    attempt = 0
    current_url = url
    response_text: str | None = None

    while attempt < MAX_RETRIES:
        attempt += 1
        # Switch to POST once we've exhausted GET attempts
        as_post = attempt >= SWITCH_TO_QLEVER_AFTER
        # Switch to QLever URL after SWITCH_TO_QLEVER_AFTER attempts
        if current_url != QLEVER_URL and attempt > SWITCH_TO_QLEVER_AFTER:
            logging.info(
                "Switching to QLever endpoint due to repeated limitations/timeouts.",
            )
            current_url = QLEVER_URL

        logging.info(
            f"Fetching query {query_file} (attempt {attempt}/{MAX_RETRIES}) via {'POST' if as_post else 'GET'} on {current_url}",
        )
        try:
            response_text = sparql_to_csv(query=query, url=current_url, as_post=as_post)
        except Exception as e:  # pragma: no cover - network/runtime errors
            logging.warning(f"Request exception: {e}")
            response_text = ""

        if _looks_timeout(response_text):
            logging.warning("Timeout detected in response text.")
        elif _looks_rate_limited(response_text):
            logging.warning("Rate limit signal detected in response text.")
        else:
            # Success path: break
            break

        # Sleep only while still below the switch threshold (original endpoint, GET phase)
        if attempt < SWITCH_TO_QLEVER_AFTER:
            logging.info(f"Retrying after {SLEEP_TIME}s...")
            sleep(SLEEP_TIME)

    if response_text is None:
        logging.error("No response received; aborting without writing file.")
        return

    if _looks_rate_limited(response_text):
        logging.error(
            "Persisting rate-limited content (all retries exhausted) – file may be incomplete.",
        )

    logging.info("Keeping unique lines (order preserving)")
    lines = response_text.splitlines()
    if not lines:
        logging.info("Empty result; nothing to write")
        return
    header = lines[0]

    # Preserve original order while removing duplicates
    seen: set[str] = set()
    unique_content_lines: list[str] = []
    for line in lines[1:]:
        if line not in seen:
            seen.add(line)
            unique_content_lines.append(line)

    processed_lines = [header]
    for line in unique_content_lines:
        processed_line = remove_wd_entity_prefix_and_Q(text=line)
        if processed_line is not None:
            processed_lines.append(processed_line)

    t_processed = "\n".join(processed_lines)

    out_path = path / output_file
    logging.info(f"Writing to {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(t_processed)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise ValueError(f"Usage: python3 {sys.argv[0]} <query_file> <output_file>")
    run(path=Path("data"), query_file=Path(sys.argv[1]), output_file=Path(sys.argv[2]))
