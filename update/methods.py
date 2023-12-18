import logging
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from update.config import MAX_WORKERS
from update.models import Task


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


def run_tasks(
    tasks: list[Task], path: Path, only=str | None, stop=str | None, skip=str | None
) -> None:
    start = time.time()

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = []
        current_group = tasks[0].group

        for task in tasks:
            if only and not task.matches_name(only):
                logging.warning(f"Skipping {task.name}")
                continue
            if stop and task.matches_name(stop):
                break
            if skip and task.matches_name(skip):
                logging.warning(f"Skipping {task.name}")
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
