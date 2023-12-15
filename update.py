#!/usr/bin/env python3
import multiprocessing
import time
from pathlib import Path

from update import (download_couples_referenced, download_dois, download_smiles, download_taxa, download_taxa_all, download_taxa_parents, download_taxa_ranks, generate_database, generate_database_chemo, generate_database_taxo)

start = time.time()

tasks = [
    ["couples_referenced", download_couples_referenced.run],
    ["dois", download_dois.run],
    ["smiles", download_smiles.run],
    ["taxa", download_taxa.run],
    ["taxa_all", download_taxa_all.run],
    ["taxa_parents", download_taxa_parents.run],
    ["taxa_ranks", download_taxa_ranks.run],
    ["generate_database_chemo", generate_database_chemo.run],
    ["generate_database_taxo", generate_database_taxo.run],
    ["generate_database", generate_database.run],
]


def run_tasks() -> None:
    for task in tasks:
        print(f"Started {task[0]}")
        start_task = time.time()
        task[1]()
        print(f" Task {task[0]} took {time.time() - start_task:.2f}s")

    print(f"Total execution time: {time.time() - start:.2f}s")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    run_tasks()
