import shutil
from pathlib import Path

from update.config import TASKS
from update.taskrunner import TaskRunner


def setup_from_fixture(tmp_path: Path):
    for file in Path("tests/fixtures").iterdir():
        shutil.copy(file, tmp_path)
    runner = TaskRunner(TASKS, tmp_path)
    runner.run_tasks(only=None, stop=None, skip="downloads", parallel=False)
