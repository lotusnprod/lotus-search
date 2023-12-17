#!/usr/bin/env python3
import argparse
import logging
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from update import download_query_as_csv, generate_database, generate_database_chemo, generate_database_taxo
from update.common import QLEVER_URL
from update.models import Group, Task

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

MAX_WORKERS = 7

start = time.time()

DownloadGroup = Group(name="downloads", parallel=True)
DatabaseGroup = Group(name="database", parallel=True)
MergingGroup = Group(name="merging", parallel=False)

TASKS = [
    Task(name="couples_referenced", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/couples_referenced.rq", "output_file": "couples.csv"}),
    Task(name="dois", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/references.rq", "output_file": "dois.csv"}),
    Task(name="smiles", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/structures.rq", "output_file": "smiles.csv"}),
    Task(name="taxa", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/taxa.rq", "output_file": "taxa.csv"}),
    Task(name="taxa_all", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/taxa_all.rq", "output_file": "taxa_all.csv", "url": QLEVER_URL}),
    Task(name="taxa_parents", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/taxa_parents.rq", "output_file": "taxa_parents.csv"}),
    Task(name="taxa_ranks", f=download_query_as_csv.run, group=DownloadGroup,
         params={"query_file": "update/queries/taxa_ranks.rq", "output_file": "taxa_ranks.csv"}),
    Task(name="generate_database_chemo", f=generate_database_chemo.run, group=DatabaseGroup),
    Task(name="generate_database_taxo", f=generate_database_taxo.run, group=DatabaseGroup),
    Task(name="generate_database", f=generate_database.run, group=MergingGroup),
]


def list_tasks(tasks: list[Task]) -> None:
    print("Tasks grouped by parallel groups:")
    group = None
    for task in tasks:
        if task.group != group:
            print("")
            group = task.group
            print(group.name, end=" > ")
        print(task.name, end=" ")
    print("")


def run_tasks(tasks: list[Task], path: Path, only=str | None, stop=str | None, skip=str | None) -> None:
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = []
        current_group = tasks[0].group

        for task in tasks:
            if only and not task.matches_name(only):
                continue
            if stop and task.matches_name(stop):
                break
            if skip and task.matches_name(skip):
                continue

            if current_group != task.group:
                for future in futures:
                    future.result()
                futures.clear()

            if current_group.parallel:
                future = executor.submit(task.run, path)
                futures.append(future)
            else:
                task.run(path)

        for future in futures:
            future.result()
    logging.info(f"Total execution time: {time.time() - start:.2f}s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", help="Only execute this task/group")
    parser.add_argument("--stop", help="Stop at this task/group")
    parser.add_argument("--skip", help="Run everything except this task/group")
    parser.add_argument("--list", action="store_true", help="List tasks")
    parser.add_argument("--data", default="data", help="Path for data, defaults to ./data")
    args = parser.parse_args()

    multiprocessing.freeze_support()
    if args.list:
        list_tasks(TASKS)
    else:
        path = Path(args.data)
        os.makedirs(path, exist_ok=True)
        run_tasks(TASKS, path=path, only=args.only, stop=args.stop, skip=args.skip)
