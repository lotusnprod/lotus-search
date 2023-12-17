from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

from update import generate_database


@patch("update.generate_database.pickle.dump")
@patch("update.generate_database.pickle.load")
@patch("update.generate_database.open", new_callable=mock_open)
def test_run_loads_and_dumps_database(mock_open, mock_pickle_load, mock_pickle_dump):
    mock_pickle_load.return_value = {"key": "value"}
    generate_database.run()
    assert mock_pickle_load.call_count == 2
    assert mock_pickle_dump.call_count == 1


@patch("update.generate_database.pickle.dump")
@patch("update.generate_database.pickle.load")
@patch("update.generate_database.open", new_callable=mock_open)
def test_run_updates_database_with_loaded_data(
    mock_open, mock_pickle_load, mock_pickle_dump
):
    mock_pickle_load.side_effect = [
        {"chemo_key": "chemo_value"},
        {"taxo_key": "taxo_value"},
    ]
    generate_database.run()
    mock_pickle_dump.assert_called_once_with(
        {"chemo_key": "chemo_value", "taxo_key": "taxo_value"},
        mock_open.return_value.__enter__.return_value,
    )
