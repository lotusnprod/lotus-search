import pytest
import requests_mock

from update.common import remove_wd_entity_prefix, wd_sparql_to_csv


def test_wd_sparql_to_csv_returns_expected_csv():
    with requests_mock.Mocker() as m:
        m.get("https://query.wikidata.org/sparql", text="expected_csv")
        result = wd_sparql_to_csv("query")
        assert result == "expected_csv"


def test_wd_sparql_to_csv_uses_provided_url():
    with requests_mock.Mocker() as m:
        m.get("https://other.url/sparql", text="expected_csv")
        result = wd_sparql_to_csv("query", "https://other.url/sparql")
        assert result == "expected_csv"


def test_remove_wd_entity_prefix_removes_prefix():
    result = remove_wd_entity_prefix("http://www.wikidata.org/entity/Q123")
    assert result == "123"


def test_remove_wd_entity_prefix_does_not_remove_other_text():
    result = remove_wd_entity_prefix("http://www.wikidata.org/entity/Q123/other")
    assert result == "123/other"
