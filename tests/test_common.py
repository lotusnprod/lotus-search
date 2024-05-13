from requests_mock import Mocker

from update.common import (
    remove_wd_entity_prefix,
    remove_wd_entity_prefix_and_Q,
    sparql_to_csv,
)


class TestWdSparqlToCsv:
    def test_returns_expected_csv(self):
        with Mocker() as m:
            m.get("https://query.wikidata.org/sparql", text="expected_csv")
            result = sparql_to_csv("query")
            assert result == "expected_csv"

    def test_returns_expected_csv_with_post(self):
        with Mocker() as m:
            m.post("https://query.wikidata.org/sparql", text="expected_csv")
            result = sparql_to_csv("query", as_post=True)
            assert result == "expected_csv"

    def test_uses_provided_url(self):
        with Mocker() as m:
            m.get("https://other.url/sparql", text="expected_csv")
            result = sparql_to_csv("query", "https://other.url/sparql")
            assert result == "expected_csv"


class TestRemoveWdEntityPrefix:
    def test_removes_prefix(self):
        result = remove_wd_entity_prefix("http://www.wikidata.org/entity/Q123")
        assert result == "Q123"

    def test_does_not_remove_other_text(self):
        result = remove_wd_entity_prefix("http://www.wikidata.org/entity/Q123/other")
        assert result == "Q123/other"

    def test_removes_prefix_2(self):
        result = remove_wd_entity_prefix_and_Q("http://www.wikidata.org/entity/Q123")
        assert result == "123"

    def test_does_not_remove_other_text_2(self):
        result = remove_wd_entity_prefix_and_Q(
            "http://www.wikidata.org/entity/Q123/other"
        )
        assert result == "123/other"

    def test_removes_prefix_3(self):
        result = remove_wd_entity_prefix_and_Q("http://www.wikidata.org/entity/P123")
        assert result == "P123"
