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

    def test_num(self):
        assert self.dm.num_taxa() == 5, f"Expected taxa, got {self.dm.num_taxa()}"
        assert (
            self.dm.num_couples() == 3
        ), f"Expected couples, got {self.dm.num_couples()}"
        assert (
            self.dm.num_structures() == 4
        ), f"Expected structures, got {self.dm.num_structures()}"

    def test_get_taxa(self):
        assert len(self.dm.get_taxa()) == 5

    def test_number_of_taxa_containing_structure(self):
        assert self.dm.get_number_of_taxa_containing_structure(1) == 2
