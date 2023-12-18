from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

from update import generate_database_chemo

#
#
# @patch("update.generate_database_chemo.ProcessPoolExecutor")
# @patch("update.generate_database_chemo.pickle.dump")
# @patch("update.generate_database_chemo.pickle.load")
# @patch("update.generate_database_chemo.open", new_callable=mock_open)
# @patch("update.generate_database_chemo.csv.reader")
# def test_run_generates_database(
#     mock_csv_reader, mock_open, mock_pickle_load, mock_pickle_dump, mock_executor
# ):
#     mock_csv_reader.return_value = iter([["c", "t", "r"]])
#     mock_pickle_load.return_value = {"key": "value"}
#     mock_executor.return_value.__enter__.return_value.map.return_value = iter(
#         [
#             (
#                 0,
#                 "smiles",
#                 "smol",
#                 "smiles_clean",
#                 "sim_fp",
#                 "sub_fp",
#                 "mol_h",
#                 "sim_fp_h",
#                 "sub_fp_h",
#             )
#         ]
#     )
#     generate_database_chemo.run()
#     assert mock_pickle_dump.call_count == 1
#
#
# @patch("update.generate_database_chemo.Chem.Mol")
# @patch("update.generate_database_chemo.Chem.MolToSmiles")
# @patch("update.generate_database_chemo.fingerprint")
# @patch("update.generate_database_chemo.Chem.PatternFingerprint")
# @patch("update.generate_database_chemo.Chem.AddHs")
# @patch("update.generate_database_chemo.standardize")
# def test_process_smiles_returns_expected_result_on_success(
#     mock_standardize,
#     mock_add_hs,
#     mock_pattern_fp,
#     mock_fingerprint,
#     mock_mol_to_smiles,
#     mock_mol,
# ):
#     mock_standardize.return_value = MagicMock()
#     mock_add_hs.return_value = MagicMock()
#     mock_mol_to_smiles.return_value = "smiles_clean"
#     mock_fingerprint.return_value = "sim_fp"
#     mock_pattern_fp.return_value = "sub_fp"
#     mock_mol.return_value.ToBinary.return_value = "mol_h"
#     result = generate_database_chemo.process_smiles((0, "smiles"))
#     assert result == (
#         0,
#         "smiles",
#         "smol",
#         "smiles_clean",
#         "sim_fp",
#         "sub_fp",
#         "mol_h",
#         "sim_fp",
#         "sub_fp",
#     )
#
#


def test_process_smiles_returns_none_on_failure():
    result = generate_database_chemo.process_smiles((0, "invalid_smiles"))
    assert result is None
