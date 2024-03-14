import pytest

from api.api import search_structures
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
class TestApiStructures:
    async def test_search_structures_pure_structure(self, data_model):
        item = Item(structure={"molecule": "C"}, limit=10, modeEnum="objects")
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[3].smiles == "C"
        assert result.description == "Structures matching the query"

    async def test_search_structures_pure_structure_ids(self, data_model):
        item = Item(structure={"molecule": "C"}, limit=10, modeEnum="ids")
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects is None

    async def test_search_error_giving_both_structure_and_formuka(self, data_model):
        item = Item(
            structure={"molecule": "C", "formula": "CH4"}, limit=10, modeEnum="objects"
        )
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)

    async def test_search_structures_formula(self, data_model):
        item = Item(structure={"formula": "CH4"}, limit=10, modeEnum="objects")
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[3].smiles == "C"
        assert result.description == "Structures matching the query"

    async def test_search_structures_formula_empty(self, data_model):
        item = Item(structure={"formula": "CH3"}, limit=10, modeEnum="ids")
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0
        assert result.description == "Structures matching the query"

    async def test_search_error_giving_both_structure_and_wid(self, data_model):
        item = Item(structure={"molecule": "C", "wid": 1})
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)

    async def test_search_invalid_structure(self, data_model):
        item = Item(structure={"molecule": "C1"}, limit=10)
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)
        item.structure.option.substructure_search = True
        with pytest.raises(Exception):
            await search_structures(item=item, dm=data_model)

    async def test_search_structures_by_substructure(self, data_model):
        item = Item(
            structure={"molecule": "C", "option": {"substructure_search": True}},
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 4
        assert result.objects[1].smiles == "C[C@H](N)O"
        assert result.objects[3].smiles == "C"

    async def test_search_structures_by_substructure_explicit_h(self, data_model):
        item = Item(
            structure={
                "molecule": "C([H])([H])([H])",
                "option": {"substructure_search": True},
            },
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 3
        assert result.objects[1].smiles == "C[C@H](N)O"
        assert result.objects[3].smiles == "C"

    async def test_search_structures_by_similarity_explicit_h(self, data_model):
        item = Item(
            structure={
                "molecule": "C([H])([H])([H])([H])",
                "option": {"substructure_search": False},
            },
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert result.objects[3].smiles == "C"

    # TODO THESE ARE NOT WORKING FOR NOW
    # async def test_search_structures_by_descriptor_simple(self, data_model):
    #     item = Item(
    #         structure={
    #             "option": {
    #                 "descriptors": {
    #                     "NumHAcceptors_min": 1,
    #                 }
    #             },
    #         },
    #         limit=10,
    #         modeEnum="objects",
    #     )
    #     result = await search_structures(item=item, dm=data_model)
    #     assert result.count == 3

    # async def test_search_structures_by_descriptor_double(self, data_model):
    #     item = Item(
    #         structure={
    #             "option": {
    #                 "descriptors": {
    #                     "NumHAcceptors_min": 1,
    #                     "NumHeteroatoms_max": 2,
    #                 }
    #             },
    #         },
    #         limit=10,
    #         modeEnum="objects",
    #     )
    #     result = await search_structures(item=item, dm=data_model)
    #     assert result.count == 2

    # async def test_search_structures_by_descriptor_triple(self, data_model):
    #     item = Item(
    #         structure={
    #             "option": {
    #                 "descriptors": {
    #                     "NumHAcceptors_min": 1,
    #                     "NumHeteroatoms_max": 2,
    #                     "MolLogP_max": -1,
    #                 }
    #             },
    #         },
    #         limit=10,
    #         modeEnum="objects",
    #     )
    #     result = await search_structures(item=item, dm=data_model)
    #     assert result.count == 0

    # async def test_search_structures_descriptors(self, data_model):
    #     item = Item(
    #         structure={
    #             "wid": 1,
    #             "option": {
    #                 "return_descriptors": True,
    #             },
    #         },
    #         limit=10,
    #         modeEnum="objects",
    #     )
    #     result = await search_structures(item=item, dm=data_model)
    #     assert result.count == 1
    #     assert result.objects[1].descriptors.MaxAbsEStateIndex == 7.833333333333334

    async def test_search_structures_sdf(self, data_model):
        item = Item(
            structure={
                "molecule": "C([H])([H])([H])([H])",
                "option": {
                    "sdf": True,
                    "substructure_search": False,
                },
            },
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert (
            result.sdf
            == "\n     RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <WID>  (3) \n3\n\n"
        )

    async def test_search_structures_sdf_multiple(self, data_model):
        item = Item(
            structure={
                "molecule": "C([H])([H])([H])",
                "option": {"sdf": True, "substructure_search": True},
            },
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert (
            result.sdf
            == "\n     RDKit          2D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5981   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1\n  2  3  1  0\n  2  4  1  0\nM  END\n>  <WID>  (1) \n1\n\n\n     RDKit          2D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5981   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  6\n  2  3  1  0\n  2  4  1  0\nM  END\n>  <WID>  (2) \n2\n\n\n     RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <WID>  (3) \n3\n\n"
        )

    async def test_search_structures_by_substructure_limits(self, data_model):
        item = Item(
            structure={"molecule": "C", "option": {"substructure_search": True}},
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 4
        item.limit = 0
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 4
        item.limit = 1
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1
        assert len(result.objects) == 1

    async def test_search_structure_restrict_taxon_exists(self, data_model):
        # We search for a compound that exist in this taxon
        item = Item(
            structure={"molecule": "CC(N)O"},
            taxon={"wid": 1},
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1

    async def test_search_structure_restrict_taxon_not_exists(self, data_model):
        # We search for a compound that does not exist in this taxon
        item = Item(structure={"molecule": "CC(N)O"}, taxon={"wid": 3}, limit=10)
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_exists(self, data_model):
        # We search for a compound that exist in this taxon
        item = Item(
            structure={"molecule": "CC(N)O"},
            reference={"wid": 1},
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1

    async def test_search_structure_restrict_reference_not_matching(self, data_model):
        # We search for a compound that does not exist in this ref
        item = Item(structure={"molecule": "CC(N)O"}, reference={"wid": 4})
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_not_existing(self, data_model):
        # We search for a compound that does not exist in this ref
        item = Item(structure={"molecule": "CC(N)O"}, reference={"wid": 666})
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_and_taxon_both_exist(
        self, data_model
    ):
        # We search for a compound that exist in this taxon
        item = Item(
            structure={"molecule": "CC(N)O"},
            reference={"wid": 1},
            taxon={"wid": 1},
            limit=10,
            modeEnum="objects",
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 1

    async def test_search_structure_restrict_reference_and_taxon_none_exist(
        self, data_model
    ):
        # We search for a compound that exist in this taxon
        item = Item(
            structure={"molecule": "CC(N)O"},
            reference={"wid": 4},
            taxon={"wid": 3},
            limit=10,
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0

    async def test_search_structure_restrict_reference_and_taxon_one_exist(
        self, data_model
    ):
        # We search for a compound that exist in this taxon
        item = Item(
            structure={"molecule": "CC(N)O"},
            reference={"wid": 1},
            taxon={"wid": 3},
            limit=10,
        )
        result = await search_structures(item=item, dm=data_model)
        assert result.count == 0
