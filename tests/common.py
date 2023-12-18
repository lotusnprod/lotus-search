import shutil
from pathlib import Path

from update.config import TASKS
from update.methods import run_tasks


def setup_from_fixture(tmp_path: Path):
    for file in Path("tests/fixtures").iterdir():
        shutil.copy(file, tmp_path)
    run_tasks(TASKS, tmp_path, only=None, stop=None, skip="downloads")


def teardown(tmp_path: Path):
    shutil.rmtree(tmp_path)
