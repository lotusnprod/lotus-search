import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict


@dataclass
class Group:
    name: str
    parallel: bool


@dataclass
class Task:
    name: str
    f: Callable | None
    group: Group
    params: Dict[str, str] = field(default_factory=dict)

    def matches_name(self, name):
        return self.name == name or self.group.name == name

    def run(self, path: Path):
        if self.f is None:
            return

        logging.info(f"Started {self.name} group={self.group.name}")
        start_task = time.time()
        self.f(**{**self.params, **{"path": path}})
        logging.info(f" Task {self.name} took {time.time() - start_task:.2f}s")
