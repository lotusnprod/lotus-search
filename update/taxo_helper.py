import csv
import logging
from collections import defaultdict, deque
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# See https://www.wikidata.org/wiki/Q2576881


def generate_taxon_parents_with_distance(path: Path) -> list[tuple[int, int, int, int]]:
    graph = defaultdict(list)

    with open(path / "taxa_parents.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                taxon_id = int(row["taxon"])
                parent_id = int(row["parent"])
                graph[taxon_id].append(parent_id)
            except (ValueError, KeyError):
                logging.exception(f"Invalid row: {row}")
                continue

    distances = []
    distance_id = 1

    # Good ol' BFS
    for source_node in graph:
        visited = {}
        queue = deque([(source_node, 0)])  # (node, distance)
        visited[source_node] = 0

        while queue:
            current_node, current_distance = queue.popleft()

            for neighbor in graph.get(current_node, []):
                if neighbor not in visited:
                    next_distance = current_distance + 1
                    visited[neighbor] = next_distance
                    queue.append((neighbor, next_distance))
                    distances.append((
                        distance_id,
                        source_node,
                        neighbor,
                        next_distance,
                    ))
                    distance_id += 1

    return distances
