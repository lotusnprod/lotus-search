import argparse
import logging
import multiprocessing
import os
from pathlib import Path

from update.config import TASKS
from update.taskrunner import TaskRunner

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", help="Only execute this task/group")
    parser.add_argument("--stop", help="Stop at this task/group")
    parser.add_argument("--skip", help="Run everything except this task/group")
    parser.add_argument("--list", action="store_true", help="List tasks")
    parser.add_argument(
        "--data",
        default="data",
        help="Path for data, defaults to ./data",
    )
    args = parser.parse_args()

    multiprocessing.freeze_support()
    runner = TaskRunner(TASKS, Path(args.data))
    if args.list:
        runner.list_tasks()
    else:
        path = Path(args.data)
        os.makedirs(path, exist_ok=True)
        runner.run_tasks(only=args.only, stop=args.stop, skip=args.skip)
