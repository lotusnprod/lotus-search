from pathlib import Path
from unittest.mock import mock_open, patch

import pytest

from update.common import remove_wd_entity_prefix, wd_sparql_to_csv
from update.download_couples_referenced import run


@patch("update.download_couples_referenced.wd_sparql_to_csv")
@patch("builtins.open", new_callable=mock_open)
def test_run_writes_to_file_on_successful_query(mock_open, mock_wd_sparql_to_csv):
    mock_wd_sparql_to_csv.return_value = "query_result"
    run(Path("."))
    mock_open.assert_called_once_with(Path(".") / "couples.csv", "w")


@patch("update.download_couples_referenced.wd_sparql_to_csv")
def test_run_retries_on_timeout(mock_wd_sparql_to_csv):
    mock_wd_sparql_to_csv.side_effect = [
        "java.util.concurrent.TimeoutException",
        "query_result",
    ]
    run(Path("."))
    assert mock_wd_sparql_to_csv.call_count == 2


@patch("update.download_couples_referenced.wd_sparql_to_csv")
def test_run_raises_exception_after_max_retries(mock_wd_sparql_to_csv):
    mock_wd_sparql_to_csv.return_value = "java.util.concurrent.TimeoutException"
    with pytest.raises(TimeoutError):
        run(Path("."), retry=0)
