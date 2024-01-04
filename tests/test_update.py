import pickle

import pytest

from model.data_model import DataModel
from tests.common import setup_from_fixture

EXPECTED_KEYS_CHEMO = [
    "structure_wid",
    "structure_sim_fps",
    "structure_sim_h_fps",
    "structure_library",
    "structure_library_h",
    "structure_id",
    "structure_ranges",
]

@pytest.fixture
def data_model(tmp_path):
    setup_from_fixture(tmp_path)
    return DataModel(tmp_path)


class TestUpdate:
    def test_pkls_exist(self, data_model):
        assert (data_model.path / "database_chemo.pkl").exists()
        assert (data_model.path / "lotus.sdf").exists()

    def test_triplets(self):
        pass
        # TODO we need to test storage entirely now

    def test_chemo(self, data_model):
        with open(data_model.path / "database_chemo.pkl", "rb") as f:
            db = pickle.load(f)
            assert len(db) == len(EXPECTED_KEYS_CHEMO)
            for expected_key in EXPECTED_KEYS_CHEMO:
                assert expected_key in db
            # We likely want to test the content as well

    def test_sdf(self, data_model):
        generated_sdf_path = data_model.path / "lotus.sdf"
        fixture_sdf_path = "tests/fixtures/lotus.sdf"

        with open(generated_sdf_path, "r") as f:
            sdf = f.read()
        # If you need to update, set r to w below, and do fo.write(sdf)
        # don't forget to remove it after
        # with open(fixture_sdf_path, "w") as fo: ## Promise I'll remove that before committing
        #    sdf_fixture = fo.write(sdf)
        with open(fixture_sdf_path, "r") as fo:
            sdf_fixture = fo.read()

        assert (
            sdf == sdf_fixture
        ), f"Content mismatch between {generated_sdf_path} and {fixture_sdf_path}"
