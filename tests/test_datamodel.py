import pytest

from model import DataModel
from tests.common import setup_from_fixture, teardown


class TestDataModel:
    tmp_path = None

    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.tmp_path = tmp_path
        setup_from_fixture(tmp_path)
        self.dm = DataModel(self.tmp_path)

    def teardown(self):
        teardown(self.tmp_path)

    def test_get_taxa(self):
        assert len(self.dm.get_taxa()) == 5
