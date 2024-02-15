import pickle

import pytest

from model.data_model import DataModel
from sdf_helpers import find_structures_bytes_ranges, mmap_file, read_selected_ranges
from tests.common import setup_from_fixture

EXPECTED_KEYS_CHEMO = [
    "structure_wid",
    "structure_sim_fps",
    "structure_sim_h_fps",
    "structure_library",
    "structure_library_h",
    "structure_id",
]


@pytest.fixture
def data_model(tmp_path):
    setup_from_fixture(tmp_path)
    return DataModel(tmp_path)


class TestUpdate:
    def test_triplets(self):
        pass
        # TODO we need to test storage entirely now

    def test_pkl_exist(self, data_model):
        assert (data_model.path / "database_chemo.pkl").exists()

    def test_pkl_headers(self, data_model):
        with open(data_model.path / "database_chemo.pkl", "rb") as f:
            db = pickle.load(f)
            assert len(db) == len(EXPECTED_KEYS_CHEMO)
            for expected_key in EXPECTED_KEYS_CHEMO:
                assert expected_key in db
            # We likely want to test the content as well

    def test_sdf_exist(self, data_model):
        assert (data_model.path / "lotus.sdf").exists()

    def test_sdf(self, data_model):
        generated_sdf_path = data_model.path / "lotus.sdf"

        # Not working anymore because of the append mechanism
        # fixture_sdf_path = "tests/fixtures/lotus.sdf"
        # with open(generated_sdf_path, "r") as f:
        #     sdf = f.read()
        # # If you need to update, set r to w below, and do fo.write(sdf)
        # # don't forget to remove it after
        # # with open(fixture_sdf_path, "w") as fo: ## Promise I'll remove that before committing
        # #    sdf_fixture = fo.write(sdf)
        # with open(fixture_sdf_path, "r") as fo:
        #     sdf_fixture = fo.read()

        # assert (
        #     sdf == sdf_fixture
        # ), f"Content mismatch between {generated_sdf_path} and {fixture_sdf_path}"

        mmaped_sdf_generated = mmap_file(generated_sdf_path)

        ranges = find_structures_bytes_ranges(mmaped_sdf_generated)
        # Hard to find bytes manually
        ranges_expected = {
            1: (0, 409),
            2: (414, 823),
            3: (828, 988),
            4: (993, 1153),
            6: (1158, 1567),
        }
        assert (
            ranges == ranges_expected
        ), f"Content mismatch between {ranges} and {ranges_expected}"

        ranges_to_read = [ranges[key] for key in list(ranges.keys())]
        block = read_selected_ranges(mmaped_sdf_generated, [ranges_to_read[2]])
        block_expected = [
            """\n     RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <WID>  (3) \n3\n\n"""
        ]
        assert (
            block == block_expected
        ), f"Content mismatch between {block} and {block_expected}"

    def test_table_exist(self, data_model):
        assert (data_model.path / "structures_table.csv").exists()

    def test_table(self, data_model):
        with open(data_model.path / "structures_table.csv", "rb") as f:
            # TODO test content
            pass

    def test_blocks_table_exist(self, data_model):
        assert (data_model.path / "structures_table.csv").exists()

    def test_blocks(self, data_model):
        with open(data_model.path / "structures_blocks_table.csv", "rb") as f:
            # TODO test content
            pass

    def test_descriptors_exist(self, data_model):
        assert (data_model.path / "descriptors_rdkit.csv").exists()

    def test_descriptors(self, data_model):
        with open(data_model.path / "descriptors_rdkit.csv", "rb") as f:
            # TODO test content
            pass

    def test_processed_exist(self, data_model):
        assert (data_model.path / "smiles_processed.csv").exists()

    def test_processed(self, data_model):
        with open(data_model.path / "smiles_processed.csv", "rb") as f:
            # TODO test content
            pass
