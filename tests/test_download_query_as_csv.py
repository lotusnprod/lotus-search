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
            with patch("update.download_query_as_csv.sleep") as mock_sleep:
                run(self.path, self.query_file, self.output_file)
                assert mock_sparql_to_csv.call_count == 2
                assert self.output_file.read_text() == "valid result"
                sleep_times = [arg_list[0][0] for arg_list in mock_sleep.call_args_list]
                assert sleep_times == [5]

    def test_retries_on_timeout_switching_endpoint(self):
        with patch("update.download_query_as_csv.sparql_to_csv") as mock_sparql_to_csv:
            mock_sparql_to_csv.side_effect = [
                "java.util.concurrent.TimeoutException",
                "java.util.concurrent.TimeoutException",
                "java.util.concurrent.TimeoutException",
                "java.util.concurrent.TimeoutException",
                "valid result",
            ]
            with patch("update.download_query_as_csv.sleep") as mock_sleep:
                run(self.path, self.query_file, self.output_file, url="foo")
                assert mock_sparql_to_csv.call_count == 5
                urls = [
                    arg_list[1]["url"] for arg_list in mock_sparql_to_csv.call_args_list
                ]
                as_post = [
                    arg_list[1]["as_post"]
                    for arg_list in mock_sparql_to_csv.call_args_list
                ]
                assert as_post == [False, False, True, True, True]
                assert (
                    urls
                    == ["foo"] * 3
                    + ["https://qlever.cs.uni-freiburg.de/api/wikidata"] * 2
                )
                assert self.output_file.read_text() == "valid result"
                assert mock_sleep.call_count == 2
                sleep_times = [arg_list[0][0] for arg_list in mock_sleep.call_args_list]
                assert sleep_times == [
                    5,
                    5,
                ]

    # TODO fix it, not working anymore
    # def test_failure_on_timeout(self):
    #     with patch("update.download_query_as_csv.sparql_to_csv") as mock_sparql_to_csv:
    #         mock_sparql_to_csv.side_effect = [
    #             "java.util.concurrent.TimeoutException",
    #             "java.util.concurrent.TimeoutException",
    #             "java.util.concurrent.TimeoutException",
    #         ]
    #         with patch("update.download_query_as_csv.sleep") as mock_sleep:
    #             with pytest.raises(TimeoutError):
    #                 run(self.path, self.query_file, self.output_file)
    #             assert mock_sparql_to_csv.call_count == 3
    #             urls = [
    #                 arg_list[1]["url"] for arg_list in mock_sparql_to_csv.call_args_list
    #             ]
    #             assert urls == ["https://query.wikidata.org/sparql"] * 3
    #             assert mock_sleep.call_count == 2
    #             sleep_times = [arg_list[0][0] for arg_list in mock_sleep.call_args_list]
    #             assert sleep_times == [5, 5]

    def test_writes_expected_result(self):
        with patch(
            "update.download_query_as_csv.sparql_to_csv"
        ) as mock_sparql_to_csv, patch(
            "update.download_query_as_csv.remove_wd_entity_prefix_and_Q"
        ) as mock_remove_wd_entity_prefix:
            mock_sparql_to_csv.return_value = "valid result"
            mock_remove_wd_entity_prefix.return_value = "expected result"
            run(self.path, self.query_file, self.output_file)
            assert self.output_file.read_text() == "expected result"
