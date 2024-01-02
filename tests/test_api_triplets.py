import pytest

from api.api import search_triplets
from api.models import Item, TripletResult
from .common import data_model


@pytest.mark.usefixtures("data_model")
class TestApiTriplets:
    @pytest.mark.asyncio
    async def test_triplets_all(self, data_model):
        item = Item()
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 6
        assert len(result.references) == 4
        assert len(result.structures) == 4
        assert len(result.taxa) == 3
        assert result.description == "Triplets matching the query"

    @pytest.mark.asyncio
    async def test_triplets_limits(self, data_model):
        item = Item(limit=0)
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 6
        item.limit = 3
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 3

    @pytest.mark.asyncio
    async def test_triplets_one_reference(self, data_model):
        item = Item(reference_wid=3)
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 2
        assert len(result.references) == 1
        assert len(result.structures) == 2
        assert len(result.taxa) == 1
        assert result.description == "Triplets matching the query"

    @pytest.mark.asyncio
    async def test_triplets_one_reference_one_compound(self, data_model):
        item = Item(reference_wid=3, structure_wid=2)
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 1
        assert len(result.references) == 1
        assert len(result.structures) == 1
        assert len(result.taxa) == 1
        assert result.description == "Triplets matching the query"

    @pytest.mark.asyncio
    async def test_triplets_one_reference_one_non_existing_compound(self, data_model):
        item = Item(reference_wid=3, structure_wid=666)
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 0
        assert len(result.references) == 0
        assert len(result.structures) == 0
        assert len(result.taxa) == 0
        assert result.description == "Triplets matching the query"

    @pytest.mark.asyncio
    async def test_triplets_one_taxon(self, data_model):
        item = Item(taxon_wid=1)
        result: TripletResult = await search_triplets(item=item, dm=data_model)
        assert result.count == 3
        assert len(result.references) == 2
        assert len(result.structures) == 2
        assert len(result.taxa) == 1
        assert result.description == "Triplets matching the query"
