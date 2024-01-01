#!/usr/bin/env python3
import csv
import logging
from collections import defaultdict, deque
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


# See https://www.wikidata.org/wiki/Q2576881
def convert_to_int_safe(s: str) -> int | None:
    try:
        result = int(s)
        return result
    except ValueError:
        logging.error(f"{s} is not a valid integer.")
        return None


def generate_taxon_parents_with_distance(path: Path, taxa: set[int]) -> list[tuple[int, int, int]]:
    graph = defaultdict(list)
    distances = []
    with open(path / "taxa_parents.csv", "r") as f:
        reader = csv.reader(f)
        headers = next(reader)
        taxon_index = headers.index("taxon")
        parent_index = headers.index("parent")

        for row in reader:
            taxon_id = convert_to_int_safe(row[taxon_index])
            parent_id = convert_to_int_safe(row[parent_index])

            if taxon_id is None or parent_id is None:
                continue
            graph[taxon_id].append(parent_id)
    # Good ol' BFS
    for node in taxa:
        visited = {node: 0}
        queue = deque([node])
        while queue:
            current_node = queue.popleft()
            current_distance = visited[current_node]

            for neighbor in graph[current_node]:
                if neighbor not in visited:
                    queue.append(neighbor)
                    visited[neighbor] = current_distance + 1
                    distances.append((node, neighbor, current_distance + 1))

    return distances
