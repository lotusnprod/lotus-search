import pytest

from api.api import search_taxa
from api.models import Item
from model.model import DataModel
from tests.common import setup_from_fixture


@pytest.fixture
def data_model(tmp_path):
    setup_from_fixture(tmp_path)
    return DataModel(tmp_path)


class TestApiTaxa:
    @pytest.mark.asyncio
    async def test_taxa_simple(self, data_model):
        item = Item(taxon_name="Taxon 1")
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1
        assert result.taxa[1].name == "Taxon 1"
        assert result.description == "Taxa matching the query"

    @pytest.mark.asyncio
    async def test_taxa_limit(self, data_model):
        expected_count = 5
        item = Item(taxon_name="Taxon")
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == expected_count
        item.limit = 0
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == expected_count
        item.limit = 1
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1

    @pytest.mark.asyncio
    async def test_taxa_restrict_structure_existing(self, data_model):
        item = Item(taxon_name="Taxon 1", structure_wid=1)
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1
        assert result.taxa[1].name == "Taxon 1"
        assert result.description == "Taxa matching the query"

    @pytest.mark.asyncio
    async def test_taxa_restrict_structure_not_existing(self, data_model):
        item = Item(taxon_name="Taxon 1", structure_wid=4)
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 0

    @pytest.mark.asyncio
    async def test_taxa_restrict_reference_existing(self, data_model):
        item = Item(taxon_name="Taxon 1", reference_wid=1)
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 1
        assert result.taxa[1].name == "Taxon 1"
        assert result.description == "Taxa matching the query"

    @pytest.mark.asyncio
    async def test_taxa_restrict_reference_not_matching(self, data_model):
        item = Item(taxon_name="Taxon 1", reference_wid=4)
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 0

    @pytest.mark.asyncio
    async def test_taxa_restrict_reference_not_existing(self, data_model):
        item = Item(taxon_name="Taxon 1", reference_wid=666)
        result = await search_taxa(item=item, dm=data_model)
        assert result.count == 0
