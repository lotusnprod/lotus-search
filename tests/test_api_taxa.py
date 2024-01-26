import pytest

from api.api import search_taxa
from api.models import (
    Filter,
    Item,
    Reference,
    StructuralFilter,
    Structure,
    TaxalFilter,
    Taxon,
)

from .common import data_model


@pytest.mark.usefixtures("data_model")
class TestApiTaxa:
    async def test_taxa_simple(self, data_model):
        item = Item(taxon={"name": "Taxon 1"})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1
        assert result.taxa[1].name == "Taxon 1"
        assert result.description == "Taxa matching the query"

    async def test_taxa_limit(self, data_model):
        expected_count = 5
        item = Item(taxon={"name": "Taxon"})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == expected_count
        item.filter.limit = 0
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == expected_count
        item.filter.limit = 1
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1

    async def test_taxa_restrict_structure_existing(self, data_model):
        item = Item(taxon={"name": "Taxon 1"}, structure={"wid": 1})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1
        assert result.taxa[1].name == "Taxon 1"
        assert result.description == "Taxa matching the query"

    async def test_taxa_restrict_structure_not_existing(self, data_model):
        item = Item(taxon={"name": "Taxon 1"}, structure={"wid": 4})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 0

    async def test_taxa_restrict_reference_existing(self, data_model):
        item = Item(taxon={"name": "Taxon 1"}, reference={"wid": 1})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1
        assert result.taxa[1].name == "Taxon 1"
        assert result.description == "Taxa matching the query"

    async def test_taxa_restrict_reference_not_matching(self, data_model):
        item = Item(taxon={"name": "Taxon 1"}, reference={"wid": 4})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 0

    async def test_taxa_restrict_reference_not_existing(self, data_model):
        item = Item(taxon={"name": "Taxon 1"}, reference={"wid": 666})
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 0
