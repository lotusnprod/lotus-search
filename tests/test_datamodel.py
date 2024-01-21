from unittest import mock

import pytest

from .common import data_model, setup_from_fixture


@pytest.mark.usefixtures("data_model")
class TestDataModel:
    def test_make_coverage_happy(self, tmp_path):
        setup_from_fixture(tmp_path)

    # Eventually TODO add taxa_names_com

    def test_get_taxon_name_from_tid(self, data_model):
        assert data_model.get_taxon_name_from_tid(1) == "Taxon 1"
        assert data_model.get_taxon_name_from_tid(666) is None

    def test_get_rank_name_from_wid(self, data_model):
        assert data_model.get_rank_name_from_wid(100) == "genus"
        assert data_model.get_rank_name_from_wid(666) is None

    def test_resolve_taxon(self, data_model):
        out = data_model.resolve_taxon("Gentiana luthea")
        assert out["names"][0]["results"][0]["currentCanonicalFull"] == "Gentiana lutea"

    def test_resolve_taxon_with_error(self, data_model):
        with mock.patch("requests.post") as mock_post:
            mock_post.side_effect = Exception("Connection error")
            with mock.patch("model.data_model.log.error") as mock_log_error:
                result = data_model.resolve_taxon("Gentiana luthea")
                mock_log_error.assert_called_once()
                assert result is None

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
        assert len(data_model.get_structures_of_reference(1)) == 2
        assert len(data_model.get_structures_of_reference(666)) == 0

    def test_get_taxa_of_structure(self, data_model):
        assert len(data_model.get_taxa_of_structure(1)) == 2
        assert len(data_model.get_taxa_of_structure(666)) == 0

    def test_get_structures_of_taxon(self, data_model):
        assert len(data_model.get_structures_of_taxon(2)) == 3
        assert len(data_model.get_structures_of_taxon(666)) == 0
        assert len(data_model.get_structures_of_taxon(2, recursive=False)) == 2

    def test_get_taxa_with_name_matching(self, data_model):
        assert len(data_model.get_taxa_with_name_matching("Taxon 1")) == 1
        assert len(data_model.get_taxa_with_name_matching("taxon 1")) == 1
        assert len(data_model.get_taxa_with_name_matching("Taxon 666")) == 0
        assert len(data_model.get_taxa_with_name_matching("Taxon")) == 5

        assert len(data_model.get_taxa_with_name_matching("Taxon 1", exact=True)) == 1
        assert len(data_model.get_taxa_with_name_matching("taxon 1", exact=True)) == 0
        assert len(data_model.get_taxa_with_name_matching("Taxon 666", exact=True)) == 0
        assert len(data_model.get_taxa_with_name_matching("Taxon", exact=True)) == 0

    def test_tsv_of_scores(self, data_model):
        assert (
            data_model.structure_get_tsv_from_scores([1, 2, 3], [0.1, 0.2, 0.3])
            == """Wikidata link	Similarity	Smiles
http://www.wikidata.org/entity/Q1	0.100	C[C@H](N)O
http://www.wikidata.org/entity/Q2	0.200	C[C@@H](N)O
http://www.wikidata.org/entity/Q3	0.300	C
"""
        )
