import requests_mock

from update.common import remove_wd_entity_prefix, sparql_to_csv


class TestWdSparqlToCsv:
    def test_returns_expected_csv(self):
        with requests_mock.Mocker() as m:
            m.get("https://query.wikidata.org/sparql", text="expected_csv")
            result = sparql_to_csv("query")
            assert result == "expected_csv"

    def test_uses_provided_url(self):
        with requests_mock.Mocker() as m:
            m.get("https://other.url/sparql", text="expected_csv")
            result = sparql_to_csv("query", "https://other.url/sparql")
            assert result == "expected_csv"


class TestRemoveWdEntityPrefix:
    def test_removes_prefix(self):
        result = remove_wd_entity_prefix("http://www.wikidata.org/entity/Q123")
        assert result == "123"

    def test_does_not_remove_other_text(self):
        result = remove_wd_entity_prefix("http://www.wikidata.org/entity/Q123/other")
        assert result == "123/other"
