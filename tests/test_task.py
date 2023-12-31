from unittest.mock import Mock, patch

from update.models import Group, Task


class TestTask:
    def test_matches_name(self):
        group = Group(name="group", parallel=False)
        task = Task("task", None, group)
        assert task.matches_name("task")
        assert task.matches_name("group")

    @patch("update.models.logging")
    def test_run_with_method(self, mock_logging, tmp_path):
        mock_function = Mock()
        task = Task("task", mock_function, Mock(name="group"))
        task.run(tmp_path)
        mock_function.assert_called_once_with(path=tmp_path)
        assert mock_logging.info.call_count == 2

    @patch("update.models.logging")
    def test_run_without_method(self, tmp_path):
        task = Task("task", None, Mock(name="group"))
        task.run(tmp_path)
