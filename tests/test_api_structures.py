import pytest

from api.api import search_structures
from api.models import Item

from .common import data_model


@pytest.mark.usefixtures("data_model")
class TestApiStructures:
    async def test_search_structures_pure_structure(self, data_model):
        item = Item(structure="C", limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.structures[3].smiles == "C"
        assert result.description == "Structures matching the query"

    async def test_search_error_giving_both_structure_and_wid(self, data_model):
        item = Item(structure="C", structure_wid="1")
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)

    async def test_search_invalid_structure(self, data_model):
        item = Item(structure="C1", limit=10)
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)
        item.substructure_search = True
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)

    async def test_search_structures_by_substructure(self, data_model):
        item = Item(structure="C", substructure_search=True, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 4
        assert result.structures[1].smiles == "C[C@H](N)O"
        assert result.structures[3].smiles == "C"
        assert result.description == "Structures matching the query"

    async def test_search_structures_by_substructure_explicit_h(self, data_model):
        item = Item(structure="C([H])([H])([H])", substructure_search=True, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 3
        assert result.structures[1].smiles == "C[C@H](N)O"
        assert result.structures[3].smiles == "C"
        assert result.description == "Structures matching the query"

    async def test_search_structures_by_similarity_explicit_h(self, data_model):
        item = Item(
            structure="C([H])([H])([H])([H])", substructure_search=False, limit=10
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.structures[3].smiles == "C"
        assert result.description == "Structures matching the query"

    async def test_search_structures_by_substructure_limits(self, data_model):
        item = Item(structure="C", substructure_search=True)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 4
        item.limit = 0
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 4
        item.limit = 1
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert len(result.structures) == 1
        assert result.description == "Structures matching the query"

    async def test_search_structure_restrict_taxon_exists(self, data_model):
        # We search for a compound that exist in this taxon
        item = Item(structure="CC(N)O", taxon_wid=1, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.structures[1].smiles == "C[C@H](N)O"
        assert result.description == "Structures matching the query"

    async def test_search_structure_restrict_taxon_not_exists(self, data_model):
        # We search for a compound that does not exist in this taxon
        item = Item(structure="CC(N)O", taxon_wid=3, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_exists(self, data_model):
        # We search for a compound that exist in this taxon
        item = Item(structure="CC(N)O", reference_wid=1, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.structures[1].smiles == "C[C@H](N)O"
        assert result.description == "Structures matching the query"

    async def test_search_structure_restrict_reference_not_matching(self, data_model):
        # We search for a compound that does not exist in this taxon
        item = Item(structure="CC(N)O", reference_wid=4)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_not_existing(self, data_model):
        # We search for a compound that does not exist in this taxon
        item = Item(structure="CC(N)O", reference_wid=666)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_and_taxon_both_exist(
        self, data_model
    ):
        # We search for a compound that exist in this taxon
        item = Item(structure="CC(N)O", reference_wid=1, taxon_wid=1, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.structures[1].smiles == "C[C@H](N)O"
        assert result.description == "Structures matching the query"

    async def test_search_structure_restrict_reference_and_taxon_none_exist(
        self, data_model
    ):
        # We search for a compound that exist in this taxon
        item = Item(structure="CC(N)O", reference_wid=4, taxon_wid=3, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_and_taxon_one_exist(
        self, data_model
    ):
        # We search for a compound that exist in this taxon
        item = Item(structure="CC(N)O", reference_wid=1, taxon_wid=3, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0
