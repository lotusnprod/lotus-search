"""Backward-compatible access point for DataModel.

Some Dash pages and legacy code import `DataModel` from `model.model`.
This thin shim simply re-exports the implementation from `data_model`.
"""
from __future__ import annotations

from .data_model import DataModel

__all__ = ["DataModel"]

