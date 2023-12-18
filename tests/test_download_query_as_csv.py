from unittest.mock import patch

import pytest

from update.download_query_as_csv import run


class TestRunQueryToCSV:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.path = tmp_path
        self.query_file = tmp_path / "query.sparql"
        self.query_file.write_text("SELECT ?item WHERE {?item wdt:P31 wd:Q5.} LIMIT 1")
        self.output_file = tmp_path / "output.csv"

    def test_retries_on_timeout(self):
        with patch("update.download_query_as_csv.sparql_to_csv") as mock_sparql_to_csv:
            mock_sparql_to_csv.side_effect = [
                "java.util.concurrent.TimeoutException",
                "valid result",
            ]
            run(self.path, self.query_file, self.output_file)
            assert mock_sparql_to_csv.call_count == 2
            assert self.output_file.read_text() == "valid result"

    def test_writes_expected_result(self):
        with patch(
            "update.download_query_as_csv.sparql_to_csv"
        ) as mock_sparql_to_csv, patch(
            "update.download_query_as_csv.remove_wd_entity_prefix"
        ) as mock_remove_wd_entity_prefix:
            mock_sparql_to_csv.return_value = "valid result"
            mock_remove_wd_entity_prefix.return_value = "expected result"
            run(self.path, self.query_file, self.output_file)
            assert self.output_file.read_text() == "expected result"
