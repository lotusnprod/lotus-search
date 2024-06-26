#!/usr/bin/env python3
import csv
import logging
import multiprocessing
import pickle  # S403
from concurrent.futures import ProcessPoolExecutor, as_completed

# from itertools import islice
from pathlib import Path

from rdkit import RDLogger
from rdkit.Chem import Mol, rdSubstructLibrary
from tqdm import tqdm

from chemistry_helpers import process_smiles
from sdf_helpers import find_structures_bytes_ranges, mmap_file, write_mols_to_sdf

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def export_descriptors_to_csv(descriptors: dict, path: Path) -> None:
    file_exists = path.exists()
    with open(path, "a") as f:
        headers = ["smiles"] + list(descriptors[next(iter(descriptors))].keys())
        csv_writer = csv.writer(f)
        if not file_exists:
            csv_writer.writerow(headers)
        for structure, properties in descriptors.items():
            row = [structure] + [properties.get(key, "") for key in headers[1:]]
            csv_writer.writerow(row)


def load_processed_smiles(path: Path) -> set:
    processed_smiles_set: set = set()
    processed_smiles_file = path / "smiles_processed.csv"
    if processed_smiles_file.exists():
        with open(processed_smiles_file, "r") as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            processed_smiles_set.update(row[1] for row in reader)
    return processed_smiles_set


def run(path: Path) -> None:
    processed_smiles_set = load_processed_smiles(path)
    smileses = []
    links = []
    with open(path / "structures.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        # for row in islice(reader, 100):
        for row in reader:
            c, smi, cano = row
            smiles = smi.strip() if smi else cano.strip()
            if smiles not in processed_smiles_set:
                smileses.append(smiles)
                links.append(int(c))

    if not smileses:
        logging.info("No new SMILES to process. Exiting.")
        return

    max_workers = multiprocessing.cpu_count()

    smis = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()

    library = rdSubstructLibrary.SubstructLibrary(smis, fps)

    mols_h = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps_h = rdSubstructLibrary.PatternHolder()

    library_h = rdSubstructLibrary.SubstructLibrary(mols_h, fps_h)

    inchis = []
    inchikeys = []
    formulas = []
    sdf_blocks = []
    p_smileses = []
    p_smols = []
    p_sim_fps = []
    p_sim_h_fps = []
    p_links = []
    smis_no_stereo = []
    inchis_no_stereo = []
    inchikeys_no_stereo = []
    # descriptors_m = {}
    descriptors_r = {}

    logging.info("Generating the chemical libraries")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_smiles, smiles): smiles
            for smiles in enumerate(smileses)
        }
        results = tuple(
            tqdm(as_completed(futures), total=len(smileses), desc="Processing SMILES")
        )
        for future in results:
            if future is not None:
                result = future.result()
                if result is not None:
                    (
                        nid,
                        smiles,
                        smol,
                        smiles_clean,
                        inchi_clean,
                        inchikey_clean,
                        formulas_clean,
                        mol_block,
                        sim_fp,
                        sub_fp,
                        # desc_mordred,
                        desc_rdkit,
                        mol_h,
                        sim_fp_h,
                        sub_fp_h,
                        smiles_no_stereo,
                        inchi_no_stereo,
                        inchikey_no_stereo,
                    ) = result

                    mols_h.AddMol(Mol(mol_h))
                    fps_h.AddFingerprint(sub_fp_h)
                    p_sim_h_fps.append(sim_fp_h)

                    smis.AddSmiles(smiles_clean)
                    inchis.append(inchi_clean)
                    inchikeys.append(inchikey_clean)
                    formulas.append(formulas_clean)
                    sdf_blocks.append((links[nid], mol_block))
                    fps.AddFingerprint(sub_fp)
                    p_sim_fps.append(sim_fp)

                    p_smols.append(smol)
                    # TODO check which one to use
                    # See #109
                    # p_smileses.append(smiles)
                    p_smileses.append(smiles_clean)
                    smis_no_stereo.append(smiles_no_stereo)
                    inchis_no_stereo.append(inchi_no_stereo)
                    inchikeys_no_stereo.append(inchikey_no_stereo)
                    # descriptors_m[smiles_clean] = desc_mordred
                    # TODO check which one to use
                    # See #109
                    descriptors_r[smiles_clean] = desc_rdkit

                    p_links.append(links[nid])
    logging.info("Finished generating the chemical libraries")

    logging.info("Generating and exporting SDF")
    write_mols_to_sdf(path / "lotus.sdf", sdf_blocks)

    logging.info("Indexing SDF")
    mmaped_sdf = mmap_file(path / "lotus.sdf")
    structures_ranges = find_structures_bytes_ranges(mmaped_sdf)

    logging.info("Generating database")
    database = {
        "structure_wid": p_links,
        "structure_sim_fps": p_sim_fps,
        "structure_sim_h_fps": p_sim_h_fps,
        "structure_library": library.Serialize(),
        "structure_library_h": library_h.Serialize(),
        "structure_id": {i[1]: i[0] for i in enumerate(p_links)},
    }
    # print(database)

    logging.info("Exporting primary table")
    smiles_file_path = path / "structures_table.csv"
    file_exists = smiles_file_path.exists()
    with open(smiles_file_path, "a") as f:  # Append mode to avoid overwriting
        csv_file = csv.writer(f)
        if not file_exists:
            # TODO check if we want the clean SMILES or not in this table
            csv_file.writerow(
                [
                    "structure",
                    "structure_smiles",
                    "structure_smiles_no_stereo",
                    "structure_inchi",
                    "structure_inchi_no_stereo",
                    "structure_inchikey",
                    "structure_inchikey_no_stereo",
                    "structure_formula",
                ]
            )
        csv_file.writerows(
            zip(
                p_links,
                p_smileses,
                smis_no_stereo,
                inchis,
                inchis_no_stereo,
                inchikeys,
                inchikeys_no_stereo,
                formulas,
            )
        )

    # TODO this is not finished.
    logging.info("Exporting blocks table")
    smiles_file_path = path / "structures_blocks_table.csv"
    file_exists = smiles_file_path.exists()
    with open(smiles_file_path, "a") as f:  # Append mode to avoid overwriting
        csv_file = csv.writer(f)
        if not file_exists:
            csv_file.writerow(
                [
                    "structure",
                    "block_range",
                ]
            )
        csv_file.writerows(
            zip(
                structures_ranges,
                structures_ranges.values(),
            )
        )

    logging.info("Exporting rdkit descriptors")
    export_descriptors_to_csv(descriptors_r, path / "descriptors_rdkit.csv")

    # logging.info("Exporting mordred descriptors")
    # export_descriptors_to_csv(descriptors_m, path / "descriptors_mordred.csv")

    logging.info("Exporting processed smiles")
    smiles_file_path = path / "smiles_processed.csv"
    file_exists = smiles_file_path.exists()
    with open(smiles_file_path, "a") as f:  # Append mode to avoid overwriting
        csv_file = csv.writer(f)
        # Write a csv with header, structure_id and structure_smiles
        if not file_exists:
            csv_file.writerow(["structure", "structure_smiles"])
        # from the two arrays p_links and p_smileses respectively
        csv_file.writerows(zip(p_links, p_smileses))

    logging.info("Exporting database")
    with open(path / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)

    logging.info("Finished exporting")


if __name__ == "__main__":
    run(Path("data"))
