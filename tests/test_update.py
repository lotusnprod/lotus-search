import pickle
import shutil
from pathlib import Path

import pytest

from tests.common import setup_from_fixture, teardown
from update.config import TASKS
from update.methods import run_tasks

EXPECTED_KEYS_CHEMO = [
    "structure_smiles",
    "structure_wid",
    "structure_sim_fps",
    "structure_sim_h_fps",
    "structure_library",
    "structure_library_h",
    "structure_id",
    "t2c",
    "c2t",
    "tc2r",
]
EXPECTED_KEYS_TAXO = [
    "taxonomy_direct_parents",
    "taxonomy_names",
    "taxonomy_ranks",
    "taxonomy_children",
    "taxonomy_parents_with_distance",
    "taxonomy_ranks_names",
]


class TestUpdate:
    tmp_path = None

    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.tmp_path = tmp_path
        setup_from_fixture(tmp_path)

    def teardown(self):
        teardown(self.tmp_path)

    def test_pkls_exist(self):
        assert (self.tmp_path / "database_chemo.pkl").exists()
        assert (self.tmp_path / "database_taxo.pkl").exists()
        assert (self.tmp_path / "database.pkl").exists()
        assert (self.tmp_path / "lotus.sdf").exists()

    def test_chemo(self):
        with open(self.tmp_path / "database_chemo.pkl", "rb") as f:
            db = pickle.load(f)
            assert len(db) == len(EXPECTED_KEYS_CHEMO)
            for expected_key in EXPECTED_KEYS_CHEMO:
                assert expected_key in db
            # We likely want to test the content as well

    def test_taxo(self):
        with open(self.tmp_path / "database_taxo.pkl", "rb") as f:
            db = pickle.load(f)
            assert len(db) == len(EXPECTED_KEYS_TAXO)
            for expected_key in EXPECTED_KEYS_TAXO:
                assert expected_key in db
                assert len(db[expected_key]) > 0, f"Empty key: {expected_key}"
            # We likely want to test the content as well

    def test_merge(self):
        with open(self.tmp_path / "database.pkl", "rb") as f:
            db = pickle.load(f)
            assert len(db) == len(EXPECTED_KEYS_TAXO + EXPECTED_KEYS_CHEMO)
            for expected_key in EXPECTED_KEYS_TAXO + EXPECTED_KEYS_CHEMO:
                assert expected_key in db
                assert len(db[expected_key]) > 0, f"Empty key: {expected_key}"
            # For this one this is probably fine as we tested in the two others

    def test_sdf(self):
        with open(self.tmp_path / "lotus.sdf", "r") as f:
            sdf = f.read()
            # If you need to update, set r to w below, and do fo.write(sdf)
            # don't forget to remove it after
            with open("tests/fixtures/lotus.sdf", "r") as fo:
                sdf_fixture = fo.read()
                assert sdf == sdf_fixture