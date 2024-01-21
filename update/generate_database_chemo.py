#!/usr/bin/env python3
import csv
import logging
import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor
# from itertools import islice
from pathlib import Path

from rdkit import RDLogger
from rdkit.Chem import Mol, rdSubstructLibrary

from chemistry_helpers import process_smiles
from sdf_helpers import mmap_file, find_structures_bytes_ranges, write_mols_to_sdf

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def load_processed_smiles(path: Path) -> set:
    processed_smiles_set = set()
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
        # for x in islice(reader, 100):
        for x in reader:
            c, smi, cano = x
            if smi == "":
                smi = cano
            smiles = smi.strip()
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
    sdf_blocks = []
    p_smileses = []
    p_smols = []
    p_sim_fps = []
    p_sim_h_fps = []
    p_links = []
    smis_no_stereo = []
    inchis_no_stereo = []

    logging.info("Generating the chemical libraries")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(process_smiles, enumerate(smileses), chunksize=1000)
        for result in results:
            if result is not None:
                (
                    nid,
                    smiles,
                    smol,
                    smiles_clean,
                    inchi_clean,
                    inchikey_clean,
                    mol_block,
                    sim_fp,
                    sub_fp,
                    # desc_mordred,
                    # desc_rdkit,
                    mol_h,
                    sim_fp_h,
                    sub_fp_h,
                    smiles_no_stereo,
                    inchi_no_stereo,
                ) = result

                mols_h.AddMol(Mol(mol_h))
                fps_h.AddFingerprint(sub_fp_h)
                p_sim_h_fps.append(sim_fp_h)

                smis.AddSmiles(smiles_clean)
                inchis.append(inchi_clean)
                inchikeys.append(inchikey_clean)
                sdf_blocks.append((links[nid], mol_block))
                fps.AddFingerprint(sub_fp)
                p_sim_fps.append(sim_fp)

                p_smols.append(smol)
                p_smileses.append(smiles)
                smis_no_stereo.append(smiles_no_stereo)
                inchis_no_stereo.append(inchi_no_stereo)

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
        "structure_ranges": structures_ranges,
    }
    # print(database)
    # TODO add BLOCKS table based on the ranges

    # TODO decide where to put InChI(Key)s

    logging.info("Exporting processed smiles")
    smiles_file_path = path / "smiles_processed.csv"
    file_exists = smiles_file_path.exists()
    with open(smiles_file_path, "a") as f:  # Append mode to avoid overwriting
        csv_file = csv.writer(f)
        # Write a csv with header, structure_id and structure_smiles
        if not file_exists:
            csv_file.writerow(["structure", "structure_smiles"])
        # from the two arrays p_links and p_smileses  respectively
        csv_file.writerows(zip(p_links, p_smileses))

    logging.info("Exporting database")
    with open(path / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)

    logging.info("Finished exporting")


if __name__ == "__main__":
    run(Path("data"))
