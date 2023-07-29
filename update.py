import multiprocessing
from pathlib import Path
import time
from update import download_couples_referenced, download_smiles, download_taxa, generate_database

start = time.time()
DATA_PATH = Path("data")

tasks = [
    ["couples_referenced", download_couples_referenced.run],
    ["smiles", download_smiles.run],
    ["taxa", download_taxa.run],
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
