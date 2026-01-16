"""Central DataModel provider for Dash pages.

Ensures a single DataModel instance is reused across all Dash UI modules
instead of each page instantiating its own (saves memory & startup time).

This does NOT change any business logic; it simply centralizes object reuse.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from model.model import DataModel


@lru_cache(maxsize=1)
def get_data_model(path: str | None = None) -> DataModel:
    """Return a cached DataModel instance.

    Parameters
    ----------
    path: optional path to data directory (only honored on first call).
    """
    data_path = Path(path) if path else Path("data")
    return DataModel(data_path)
