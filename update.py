#!/usr/bin/env python3
import multiprocessing
import time
from pathlib import Path

from update import download_couples_referenced, download_smiles, download_taxonomy_parenting, generate_database

start = time.time()
DATA_PATH = Path("data")

tasks = [
    ["couples_referenced", download_couples_referenced.run],
    ["smiles", download_smiles.run],
    ["taxonomy_parenting", download_taxonomy_parenting.run],
    ["generate_database", generate_database.run]
]


def run_tasks() -> None:
    for task in tasks:
        print(f"Started {task[0]}")
        start_task = time.time()
        task[1](DATA_PATH)
        print(f" Task {task[0]} took {time.time() - start_task:.2f}s")

    print(f"Total execution time: {time.time() - start:.2f}s")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    run_tasks()
