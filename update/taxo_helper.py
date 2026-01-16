import csv
import logging
import os
import time
from collections import defaultdict, deque
from pathlib import Path
from typing import DefaultDict, Deque, Dict, List, Tuple, Set

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# See https://www.wikidata.org/wiki/Q2576881

__all__ = ["generate_taxon_parents_with_distance"]

# Type alias for readability
DistanceTuple = Tuple[int, int, int, int]

# Environment knobs (non-breaking defaults)
_FAST_FLAG = "LOTUS_FAST_TAXO"
_MAX_DEPTH_ENV = "LOTUS_TAXO_MAX_DEPTH"  # optional int
_SAMPLE_SKIPS_ENV = "LOTUS_TAXO_SAMPLE_SKIPS"  # number of skipped rows to log
_PROGRESS_FLAG = "LOTUS_TAXO_PROGRESS"  # progress disabled by default now
_PROGRESS_EVERY_ENV = "LOTUS_TAXO_PROGRESS_EVERY"  # sources interval (default 2000)


def _parse_max_depth() -> int | None:
    val = os.getenv(_MAX_DEPTH_ENV)
    if not val:
        return None
    try:
        d = int(val)
        return d if d >= 1 else None
    except ValueError:
        logging.warning(
            "Invalid %s=%s (must be positive int), ignoring",
            _MAX_DEPTH_ENV,
            val,
        )
        return None


def _now() -> float:
    return time.perf_counter()


def _progress_enabled() -> bool:
    # Force disabled (calculation was fast; user requested removal)
    return False


def _progress_interval() -> int:
    try:
        return max(1, int(os.getenv(_PROGRESS_EVERY_ENV, "2000")))
    except ValueError:
        return 2000


def _format_eta(elapsed: float, done: int, total: int) -> str:
    if done == 0 or total <= 0:
        return "ETA:?"
    rate = elapsed / done
    remaining = (total - done) * rate
    if remaining < 60:
        return f"ETA:{remaining:0.1f}s"
    if remaining < 3600:
        return f"ETA:{remaining / 60:0.1f}m"
    return f"ETA:{remaining / 3600:0.1f}h"


def generate_taxon_parents_with_distance(path: Path) -> List[DistanceTuple]:
    """Compute (id, source_taxon, parent_taxon, distance) tuples for all ancestry paths.

    Skips rows whose taxon/parent fields cannot be parsed as plain integers
    (e.g. internal blank node IRIs like `.well-known/genid/...`).
    Avoids raising stack traces; logs a concise summary instead.
    """
    graph: DefaultDict[int, List[int]] = defaultdict(list)

    max_depth = _parse_max_depth()
    sample_limit = 0
    try:
        sample_limit = int(os.getenv(_SAMPLE_SKIPS_ENV, "0"))
    except ValueError:
        sample_limit = 0
    skipped_samples: List[str] = []

    t_start_read = _now()
    csv_path = path / "taxa_parents.csv"
    if not csv_path.exists():
        logging.warning(
            "taxa_parents.csv not found at %s; returning empty list",
            csv_path,
        )
        return []

    processed_rows = 0
    skipped_rows = 0
    duplicate_edges = 0

    # Build adjacency list (child -> list[parent]) ensuring uniqueness per child
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        per_child_seen: Dict[int, Set[int]] = defaultdict(set)
        for row in reader:
            processed_rows += 1
            taxon_raw = row.get("taxon")
            parent_raw = row.get("parent")
            if (
                not taxon_raw
                or not parent_raw
                or not taxon_raw.isdigit()
                or not parent_raw.isdigit()
            ):
                skipped_rows += 1
                if sample_limit and len(skipped_samples) < sample_limit:
                    skipped_samples.append(str(row))
                continue
            taxon_id = int(taxon_raw)
            parent_id = int(parent_raw)
            # Skip self loops just in case
            if taxon_id == parent_id:
                skipped_rows += 1
                if sample_limit and len(skipped_samples) < sample_limit:
                    skipped_samples.append(str(row))
                continue
            if parent_id in per_child_seen[taxon_id]:
                duplicate_edges += 1
                continue
            per_child_seen[taxon_id].add(parent_id)
            graph[taxon_id].append(parent_id)

    t_after_read = _now()
    if skipped_rows:
        logging.info(
            "taxa_parents: %d rows processed, %d skipped (non-integer/invalid), %d duplicate edges ignored",
            processed_rows,
            skipped_rows,
            duplicate_edges,
        )
    if skipped_samples:
        logging.info(
            "Sample skipped rows (first %d): %s",
            len(skipped_samples),
            skipped_samples,
        )

    # FAST MODE
    if os.getenv(_FAST_FLAG) == "1":
        logging.info(
            "Using fast taxonomy ancestry computation (%s=1)%s",
            _FAST_FLAG,
            f", max_depth={max_depth}" if max_depth else "",
        )
        t_fast_start = _now()
        cache: Dict[int, Dict[int, int]] = {}
        visiting: Set[int] = set()
        total_sources = len(graph)
        progress_on = _progress_enabled()
        progress_interval = _progress_interval()
        processed_sources = 0
        total_ancestors = 0
        start_time = time.perf_counter()

        def ancestors(node: int) -> Dict[int, int]:
            if node in cache:
                return cache[node]
            if node in visiting:
                logging.warning(
                    "Cycle detected at taxon %s; skipping recursive expansion",
                    node,
                )
                return {}
            visiting.add(node)
            acc: Dict[int, int] = {}
            for parent in graph.get(node, ()):  # immediate parents
                if max_depth is not None and 1 > max_depth:
                    continue
                if parent not in acc or 1 < acc[parent]:
                    acc[parent] = 1
                if max_depth is None or max_depth > 1:
                    for and, dist in ancestors(parent).items():
                        and = dist + 1
                        if max_depth is not None and and > max_depth:
                            continue
                        if and not in acc or and < acc[and]:
                            acc[and] = and
            visiting.remove(node)
            cache[node] = acc
            return acc

        distances: List[DistanceTuple] = []
        distance_id = 1
        for source_node in graph.keys():
            and_map = ancestors(source_node)
            total_ancestors += len(and_map)
            for and, dist in sorted(and_map.items(), key=lambda x: (x[1], x[0])):
                distances.append((distance_id, source_node, and, dist))
                distance_id += 1
            processed_sources += 1
            if progress_on and processed_sources % progress_interval == 0:
                elapsed = time.perf_counter() - start_time
                pct = (processed_sources / total_sources) * 100 if total_sources else 0
                eta = _format_eta(elapsed, processed_sources, total_sources)
                logging.info(
                    "TAXO_PROGRESS mode=fast sources=%d/%d (%.2f%%) ancestors=%d elapsed=%.1fs %s",
                    processed_sources,
                    total_sources,
                    pct,
                    total_ancestors,
                    elapsed,
                    eta,
                )
        t_fast_end = _now()
        if progress_on:
            elapsed = time.perf_counter() - start_time
            logging.info(
                "TAXO_PROGRESS mode=fast COMPLETE sources=%d ancestors=%d total=%.1fs",
                processed_sources,
                total_ancestors,
                elapsed,
            )
        logging.debug(
            "TAXO_TIMING fast_mode read=%.3fs compute=%.3fs total=%.3fs",
            (t_after_read - t_start_read),
            (t_fast_end - t_fast_start),
            (t_fast_end - t_start_read),
        )
        return distances

    # Original BFS implementation
    t_bfs_start = _now()
    distances: List[DistanceTuple] = []
    distance_id = 1
    total_sources = len(graph)
    progress_on = _progress_enabled()
    progress_interval = _progress_interval()
    processed_sources = 0
    total_ancestors = 0
    start_time = time.perf_counter()
    # BFS per source taxon
    for (
        source_node,
        parents,
    ) in graph.items():  # preserves insertion order of source nodes
        if not parents:
            processed_sources += 1
            continue
        prev_len = len(distances)
        visited: Dict[int, int] = {source_node: 0}
        queue: Deque[Tuple[int, int]] = deque([(source_node, 0)])
        append_distance = distances.append
        get_neighbors = graph.get
        while queue:
            current_node, current_dist = queue.popleft()
            for neighbor in get_neighbors(current_node, ()):  # parents list
                if neighbor not in visited:
                    next_distance = current_dist + 1
                    if max_depth is not None and next_distance > max_depth:
                        continue
                    visited[neighbor] = next_distance
                    queue.append((neighbor, next_distance))
                    append_distance((distance_id, source_node, neighbor, next_distance))
                    distance_id += 1
        added = len(distances) - prev_len
        total_ancestors += added
        processed_sources += 1
        if progress_on and processed_sources % progress_interval == 0:
            elapsed = time.perf_counter() - start_time
            pct = (processed_sources / total_sources) * 100 if total_sources else 0
            eta = _format_eta(elapsed, processed_sources, total_sources)
            logging.info(
                "TAXO_PROGRESS mode=bfs sources=%d/%d (%.2f%%) ancestors=%d elapsed=%.1fs %s",
                processed_sources,
                total_sources,
                pct,
                total_ancestors,
                elapsed,
                eta,
            )
    t_bfs_end = _now()
    if progress_on:
        elapsed = time.perf_counter() - start_time
        logging.info(
            "TAXO_PROGRESS mode=bfs COMPLETE sources=%d ancestors=%d total=%.1fs",
            processed_sources,
            total_ancestors,
            elapsed,
        )
    logging.debug(
        "TAXO_TIMING bfs_mode read=%.3fs compute=%.3fs total=%.3fs%s",
        (t_after_read - t_start_read),
        (t_bfs_end - t_bfs_start),
        (t_bfs_end - t_start_read),
        f", max_depth={max_depth}" if max_depth else "",
    )
    return distances
