import pytest

from api.api import search_references
from api.models import (  # ReferenceOption,
    Item,
    ReferenceItem,
    StructureItem,
    StructureOption,
    TaxonItem,
    TaxonOption,
)

from .common import data_model


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

    async def test_search_references_pure_reference(self, data_model):
        item = Item(reference={"doi": "42.1/1"}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.description == "References matching the query"

    async def test_search_references_error_giving_doi_and_id(self, data_model):
        item = Item(reference={"doi": "42.1/1", "wid": 1})
        with pytest.raises(Exception):
            await search_references(item=item, dm=data_model)

    async def test_search_references_from_taxon_structure(self, data_model):
        item = Item(taxon={"wid": "1"}, structure={"wid": 1})
        result = await search_references(item=item, dm=data_model)
        assert result.count == 2

    async def test_search_references_with_taxon_existing(self, data_model):
        item = Item(reference={"doi": "42.1/1"}, taxon={"wid": 1}, modeEnum="objects")
        result = await search_references(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[1].doi == "42.1/1"
        assert result.description == "References matching the query"

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
        assert result.objects[1].doi == "42.1/1"
        assert result.description == "References matching the query"

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
