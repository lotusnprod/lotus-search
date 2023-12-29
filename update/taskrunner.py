import logging
import time
from concurrent.futures import Future, ProcessPoolExecutor
from pathlib import Path

from update.config import MAX_WORKERS
from update.models import Task


class TaskRunner:
    def __init__(self, tasks: list[Task], path: Path):
        self.tasks = tasks
        self.path = path
        self.workers = MAX_WORKERS

    def list_tasks(self) -> None:
        print("Tasks grouped by parallel groups:")
        group = None
        for task in self.tasks:
            if task.group != group:
                print("")
                group = task.group
                print(group.name, end=" > ")
            print(task.name, end=" ")
        print("")

    def run_tasks(
        self,
        only: str | None = None,
        stop: str | None = None,
        skip: str | None = None,
        parallel: bool = True,
    ) -> None:
        """
        Run tasks in order, optionally skipping or stopping at a task.

        :param only:
        :param stop:
        :param skip:
        :param parallel: Set to False in tests so we can get accurate coverage
        :return:
        """
        start = time.time()

        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures: list[Future[None]] = []
            current_group = self.tasks[0].group
            for task in self.tasks:
                if only and not task.matches_name(only):
                    logging.warning(f"Skipping {task.name}")
                    continue
                if stop and task.matches_name(stop):
                    break
                if skip and task.matches_name(skip):
                    logging.warning(f"Skipping {task.name}")
                    continue

                if current_group != task.group:
                    for f in futures:
                        f.result()
                    futures.clear()

                if current_group.parallel and parallel:
                    future: Future[None] = executor.submit(task.run, self.path)
                    futures.append(future)
                else:
                    task.run(self.path)

            for f in futures:
                f.result()
        logging.info(f"Total execution time: {time.time() - start:.2f}s")
