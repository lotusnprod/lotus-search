import sys
from io import StringIO
from pathlib import Path
from unittest.mock import Mock, call

from update.models import Group, Task
from update.taskrunner import TaskRunner


def f1(foo: "", path: Path):  # pragma: no cover
    with open(path / "file1", "w") as f:
        f.write(foo)


def f2(foo: "", path: Path):  # pragma: no cover
    with open(path / "file2", "w") as f:
        f.write(foo)


class TestTaskRunner:
    def test_one_task(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        tasks = [
            Task(
                name="task1",
                f=mock_call.run,
                params={"foo": "bar"},
                group=Group(name="group1", parallel=False),
            ),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(parallel=True)
        assert mock_call.run.call_count == 1
        mock_call.run.assert_called_once_with(foo="bar", path=tmp_path)

    def test_multiple_tasks_same_group(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group = Group(name="group1", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(parallel=True)
        assert mock_call.run.call_count == 2
        mock_call.run.assert_has_calls([call(foo="bar", path=tmp_path), call(foo="bim", path=tmp_path)])

    def test_multiple_tasks_different_groups(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(parallel=True)
        assert mock_call.run.call_count == 2
        mock_call.run.assert_has_calls([call(foo="bar", path=tmp_path), call(foo="bim", path=tmp_path)])

    def test_multiple_tasks_only_one_task(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(only="task2", parallel=True)
        assert mock_call.run.call_count == 1
        mock_call.run.assert_called_once_with(foo="bim", path=tmp_path)

    def test_multiple_tasks_only_one_group(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(only="group2", parallel=True)
        assert mock_call.run.call_count == 1
        mock_call.run.assert_called_once_with(foo="bim", path=tmp_path)

    def test_multiple_tasks_skip_one_task(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(skip="task1", parallel=True)
        assert mock_call.run.call_count == 1
        mock_call.run.assert_called_once_with(foo="bim", path=tmp_path)

    def test_multiple_tasks_stop_at_one_task(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(stop="task2", parallel=True)
        assert mock_call.run.call_count == 1
        mock_call.run.assert_called_once_with(foo="bar", path=tmp_path)

    def test_multiple_tasks_skip_one_group(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(skip="group1", parallel=True)
        assert mock_call.run.call_count == 1
        mock_call.run.assert_called_once_with(foo="bim", path=tmp_path)

    def test_list_tasks_different_group(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        group2 = Group(name="group2", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        captured_output = StringIO()
        sys.stdout = captured_output
        runner.list_tasks()
        assert captured_output.getvalue() == "\n".join([
            "Tasks grouped by parallel groups:",
            "",
            "group1 > task1 ",
            "group2 > task2 ",
            "",
        ])

    def test_list_tasks_same_group(self, tmp_path):
        mock_call = Mock()
        mock_call.run.return_value = True
        group1 = Group(name="group1", parallel=False)
        tasks = [
            Task(name="task1", f=mock_call.run, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=mock_call.run, params={"foo": "bim"}, group=group1),
        ]
        runner = TaskRunner(tasks, tmp_path)
        captured_output = StringIO()
        sys.stdout = captured_output
        runner.list_tasks()
        assert captured_output.getvalue() == "\n".join([
            "Tasks grouped by parallel groups:",
            "",
            "group1 > task1 task2 ",
            "",
        ])

    def test_parallelism(self, tmp_path):
        # We need to have side effects here to test parallelism as it creates new processes
        group1 = Group(name="group1", parallel=True)
        group2 = Group(name="group2", parallel=True)
        tasks = [
            Task(name="task1", f=f1, params={"foo": "bar"}, group=group1),
            Task(name="task2", f=f2, params={"foo": "bim"}, group=group2),
        ]
        runner = TaskRunner(tasks, tmp_path)
        runner.run_tasks(parallel=True)
        assert Path(tmp_path / "file1").read_text() == "bar"
        assert Path(tmp_path / "file2").read_text() == "bim"
