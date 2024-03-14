import pytest

from api.api import search_references
from api.models import (
    Item,
    ReferenceItem,
    ReferenceOption,
    StructureItem,
    StructureOption,
    TaxonItem,
    TaxonOption,
)

from tests.common import data_model


@pytest.mark.usefixtures("data_model")
class TestApiReferences:
    async def test_search_with_limits(self, data_model):
        expected_total = 4
        item = Item(reference={"doi": "42"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == expected_total
        item.limit = 0
        result = await search_references(item=item, dm=data_model)
        assert result.count == expected_total
        item.limit = 1
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1

    async def test_search_with_limits_id(self, data_model):
        expected_total = 4
        item = Item(reference={"doi": "42"}, modeEnum="ids")
        result = await search_references(item=item, dm=data_model)
        assert result.count == expected_total
        assert result.objects is None

    async def test_search_references_doi(self, data_model):
        item = Item(reference={"doi": "42.1/1"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.objects[1].title == "TITLE A with Gentiana lutea"
        assert result.objects[1].date == "2010-01-01T00:00:00Z"
        assert result.objects[1].journal == "journal A"
        assert result.description == "References matching the query"

    async def test_search_references_title(self, data_model):
        item = Item(reference={"title": "TITLE A"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.objects[1].title == "TITLE A with Gentiana lutea"
        assert result.objects[1].date == "2010-01-01T00:00:00Z"
        assert result.objects[1].journal == "journal A"
        assert result.description == "References matching the query"

    async def test_search_references_date_min(self, data_model):
        item = Item(
            reference={"option": {"date_min": "2000"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 3

    async def test_search_references_date_max(self, data_model):
        item = Item(
            reference={"option": {"date_max": "2019"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 3

    async def test_search_references_date_min_ok(self, data_model):
        item = Item(
            reference={"title": "TITLE A", "option": {"date_min": "2000"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.objects[1].title == "TITLE A with Gentiana lutea"
        assert result.objects[1].date == "2010-01-01T00:00:00Z"
        assert result.objects[1].journal == "journal A"
        assert result.description == "References matching the query"

    async def test_search_references_date_min_no(self, data_model):
        item = Item(
            reference={"title": "TITLE A", "option": {"date_min": "2011"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_date_max_ok(self, data_model):
        item = Item(
            reference={"title": "TITLE A", "option": {"date_max": "2011"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.objects[1].title == "TITLE A with Gentiana lutea"
        assert result.objects[1].date == "2010-01-01T00:00:00Z"
        assert result.objects[1].journal == "journal A"
        assert result.description == "References matching the query"

    async def test_search_references_date_max_no(self, data_model):
        item = Item(
            reference={"title": "TITLE A", "option": {"date_max": "2000"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_date_min_max_ok(self, data_model):
        item = Item(
            reference={
                "title": "TITLE A",
                "option": {"date_min": "2000", "date_max": "2011"},
            },
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.objects[1].title == "TITLE A with Gentiana lutea"
        assert result.objects[1].date == "2010-01-01T00:00:00Z"
        assert result.objects[1].journal == "journal A"
        assert result.description == "References matching the query"

    async def test_search_references_date_min_max_no(self, data_model):
        item = Item(
            reference={
                "title": "TITLE A",
                "option": {"date_min": "2011", "date_max": "2024"},
            },
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_journal_all(self, data_model):
        item = Item(
            reference={"option": {"journal": "journal"}},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 4

    async def test_search_references_journal_a(self, data_model):
        item = Item(
            reference={
                "option": {"journal": "journal A"},
            },
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 2

    async def test_search_references_journal_ok(self, data_model):
        item = Item(
            reference={
                "title": "TITLE A",
                "option": {
                    "journal": "journal A",
                },
            },
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.objects[1].title == "TITLE A with Gentiana lutea"
        assert result.objects[1].date == "2010-01-01T00:00:00Z"
        assert result.objects[1].journal == "journal A"
        assert result.description == "References matching the query"

    async def test_search_references_journal_no(self, data_model):
        item = Item(
            reference={
                "title": "TITLE A",
                "option": {"journal": "journal B"},
            },
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_title_no_hit(self, data_model):
        item = Item(reference={"title": "Foo Bar"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_title_unspecific(self, data_model):
        item = Item(reference={"title": "with"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 4

    async def test_search_references_title_lactone(self, data_model):
        item = Item(reference={"title": "lactone"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 2

    async def test_search_references_error_giving_doi_and_id(self, data_model):
        item = Item(reference={"doi": "42.1/1", "wid": 1})
        with pytest.raises(Exception):
            await search_references(item=item, dm=data_model)

    async def test_search_references_error_giving_doi_and_title(self, data_model):
        item = Item(reference={"doi": "42.1/1", "title": "Foo"})
        with pytest.raises(Exception):
            await search_references(item=item, dm=data_model)

    async def test_search_references_error_giving_id_and_title(self, data_model):
        item = Item(reference={"wid": 1, "title": "Bar"})
        with pytest.raises(Exception):
            await search_references(item=item, dm=data_model)

    async def test_search_references_from_taxon_structure(self, data_model):
        item = Item(taxon={"wid": "1"}, structure={"wid": 1})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 2

    async def test_search_references_from_taxon_structure_title(self, data_model):
        item = Item(
            reference={"title": "Foo Bar"}, taxon={"wid": "1"}, structure={"wid": 1}
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_with_taxon_existing(self, data_model):
        item = Item(reference={"doi": "42.1/1"}, taxon={"wid": 1}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1

    async def test_search_references_with_taxon_non_existing(self, data_model):
        item = Item(reference={"doi": "42.1/1"}, taxon={"wid": 666})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_with_taxon_non_matching(self, data_model):
        item = Item(reference={"doi": "42.1/1"}, taxon={"wid": 4})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_with_taxon_structure_existing(self, data_model):
        item = Item(
            reference={"doi": "42.1/1"},
            taxon={"wid": 1},
            structure={"wid": 1},
            modeEnum="objects",
        )
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1

    async def test_search_references_with_taxon_non_matching_structure_existing(
        self, data_model
    ):
        item = Item(reference={"doi": "42.1/1"}, taxon={"wid": 2}, structure={"wid": 1})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_references_with_structure_wid(self, data_model):
        item = Item(structure={"wid": 1})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 3

    async def test_search_references_with_taxon_wid(self, data_model):
        item = Item(taxon={"wid": 1})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 2
