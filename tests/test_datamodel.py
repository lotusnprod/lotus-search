import pytest

from model.model import DataModel
from tests.common import setup_from_fixture


@pytest.fixture
def data_model(tmp_path):
    setup_from_fixture(tmp_path)
    return DataModel(tmp_path)


class TestDataModel:
    def test_make_coverage_happy(self, tmp_path):
        setup_from_fixture(tmp_path)

    def test_get_taxa(self, data_model):
        assert len(data_model.get_taxa()) == 5

    def test_get_taxon_name_from_tid(self, data_model):
        assert data_model.get_taxon_name_from_tid(1) == "Taxon 1"
        assert data_model.get_taxon_name_from_tid(666) is None

    def test_get_rank_name_from_wid(self, data_model):
        assert data_model.get_rank_name_from_wid(100) == "genus"
        assert data_model.get_rank_name_from_wid(666) is None

    def test_resolve_taxon(self, data_model):
        out = data_model.resolve_taxon("Gentiana luthea")
        assert out["names"][0]["results"][0]["currentCanonicalFull"] == "Gentiana lutea"

    def test_get_ranks_string(self, data_model):
        assert data_model.get_ranks_string(1) == " (genus)"
        assert data_model.get_ranks_string(4) == " (species)"
        assert data_model.get_ranks_string(666) == ""

    def test_get_structure_smiles_from_sid(self, data_model):
        assert data_model.get_structure_smiles_from_sid(1) == "C[C@H](N)O"
        assert data_model.get_structure_smiles_from_sid(666) is None

    def test_get_reference_doi_from_rid(self, data_model):
        assert data_model.get_reference_doi_from_rid(1) == "42.1/1"
        assert data_model.get_reference_doi_from_rid(666) is None

    def test_get_references_of_taxon(self, data_model):
        assert len(data_model.get_references_of_taxon(1)) == 2
        assert len(data_model.get_references_of_taxon(666)) == 0

    def test_get_references_of_structure(self, data_model):
        assert len(data_model.get_references_of_structure(1)) == 3
        assert len(data_model.get_references_of_structure(666)) == 0

    def test_get_taxa_of_reference(self, data_model):
        assert len(data_model.get_taxa_of_reference(1)) == 1
        assert len(data_model.get_taxa_of_reference(666)) == 0

    def test_get_structures_of_reference(self, data_model):
        assert len(data_model.get_structures_of_reference(1)) == 1
        assert len(data_model.get_structures_of_reference(666)) == 0

    def test_get_taxa_of_structure(self, data_model):
        assert len(data_model.get_taxa_of_structure(1)) == 2
        assert len(data_model.get_taxa_of_structure(666)) == 0

    def test_get_structures_of_taxon(self, data_model):
        assert len(data_model.get_structures_of_taxon(2)) == 3
        assert len(data_model.get_structures_of_taxon(666)) == 0
        assert len(data_model.get_structures_of_taxon(2, recursive=False)) == 2

    def test_get_taxonomic_tree(self, data_model):
        assert data_model.get_taxonomic_tree(1) == [(5, 1), (8, 2)]
        assert data_model.get_taxonomic_tree(666) == []

    def test_get_taxa_with_name_exact(self, data_model):
        assert len(data_model.get_taxa_with_name_exact("Taxon 1")) == 1
        assert len(data_model.get_taxa_with_name_exact("taxon 1")) == 0
        assert len(data_model.get_taxa_with_name_exact("Taxon 666")) == 0
        assert len(data_model.get_taxa_with_name_exact("Taxon")) == 0
