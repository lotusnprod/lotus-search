import shutil
from pathlib import Path

import pytest

from model.data_model import DataModel
from update.config import TASKS
from update.taskrunner import TaskRunner


def setup_from_fixture(tmp_path: Path):
    for file in Path("tests/fixtures").iterdir():
        shutil.copy(file, tmp_path)
    runner = TaskRunner(TASKS, tmp_path)
    runner.run_tasks(only=None, stop=None, skip="downloads", parallel=False)


@pytest.fixture(scope="session")
def data_model(tmp_path_factory):
    temp_dir = tmp_path_factory.mktemp("data_model")
    setup_from_fixture(temp_dir)
    return DataModel(temp_dir)
