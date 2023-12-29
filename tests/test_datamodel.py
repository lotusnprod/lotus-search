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

    def test_get_rank_name_from_wid(self, data_model):
        assert data_model.get_rank_name_from_wid(100) == "genus"

    def test_resolve_taxon(self, data_model):
        out = data_model.resolve_taxon("Gentiana luthea")
        assert out["names"][0]["results"][0]["currentCanonicalFull"] == "Gentiana lutea"

    def test_get_ranks_string(self, data_model):
        pass
        # TODO: This is broken and we need that for the UI
        # assert data_model.get_ranks_string(1) == "genus"
