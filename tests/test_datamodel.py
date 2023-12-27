import pytest

from model import DataModel
from tests.common import setup_from_fixture


@pytest.fixture
def data_model(tmp_path):
    setup_from_fixture(tmp_path)
    return DataModel(tmp_path)


class TestDataModel:
    def test_get_taxa(self, data_model):
        assert len(data_model.get_taxa()) == 5
