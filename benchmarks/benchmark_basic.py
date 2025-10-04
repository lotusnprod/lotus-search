"""Minimal performance benchmark harness (no functional logic changes).

Usage:
  uv run python -m benchmarks.benchmark_basic

Environment variables:
  BENCH_ITER (int): number of iterations for repeated lookups (default 200)

Outputs JSON lines (one per benchmark) for easy diffing in CI.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

from model.data_model import DataModel


def time_block(label: str, func, *args, **kwargs):  # noqa: ANN001, D401
    start = time.perf_counter()
    result = func(*args, **kwargs)
    duration = time.perf_counter() - start
    print(
        json.dumps({
            "benchmark": label,
            "seconds": duration,
            "size": getattr(result, "__len__", lambda: None)(),
        })
    )
    return result


def main() -> None:  # pragma: no cover - benchmark utility
    data_path = Path("data")
    if not data_path.exists():
        print(json.dumps({"error": "data directory missing; run update tasks first"}))
        return
    dm = DataModel(data_path)

    iterations = int(os.environ.get("BENCH_ITER", "200"))

    # Warm-up caches
    dm.get_taxa_with_name_matching("Taxon")
    dm.get_rank_name_from_wid(1)

    time_block("taxa_name_matching_first", dm.get_taxa_with_name_matching, "Taxon")
    for _ in range(iterations):
        dm.get_taxa_with_name_matching("Taxon")
    time_block("taxa_name_matching_cached", dm.get_taxa_with_name_matching, "Taxon")

    # Structure similarity sample (choose simplest known SMILES if present)
    smiles = "C"
    time_block("structure_similarity", dm.structure_search, smiles)
    time_block("structure_substructure", dm.structure_search_substructure, smiles)

    # SDF extraction for a few IDs (if present)
    existing_ids = list(dm.structures_set())[:3]
    if existing_ids:
        time_block(
            "sdf_extraction_three", dm.get_structure_sdf_from_dict_of_sids, existing_ids
        )


if __name__ == "__main__":  # pragma: no cover
    main()
