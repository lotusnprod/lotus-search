#!/usr/bin/env python3
import csv
import logging
import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import RDLogger
from rdkit.Chem import Mol, rdSubstructLibrary

from chemistry_helpers import (fingerprint, process_smiles,
                               process_smol_and_wid, standardize,
                               write_mols_to_sdf)

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    smileses = []
    links = []
    with open(path / "smiles.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            c, smi, cano = x
            if smi == "":
                smi = cano
            smileses.append(smi)
            links.append(int(c))

    max_workers = multiprocessing.cpu_count()

    smis = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()

    library = rdSubstructLibrary.SubstructLibrary(smis, fps)

    mols_h = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps_h = rdSubstructLibrary.PatternHolder()

    library_h = rdSubstructLibrary.SubstructLibrary(mols_h, fps_h)

    p_smileses = []
    p_smols = []
    p_sim_fps = []
    p_sim_h_fps = []
    p_links = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(process_smiles, enumerate(smileses), chunksize=1000)
        for result in results:
            if result is not None:
                (
                    mid,
                    smiles,
                    smol,
                    smiles_clean,
                    sim_fp,
                    sub_fp,
                    mol_h,
                    sim_fp_h,
                    sub_fp_h,
                ) = result

                mols_h.AddMol(Mol(mol_h))
                fps_h.AddFingerprint(sub_fp_h)
                p_sim_h_fps.append(sim_fp_h)

                smis.AddSmiles(smiles_clean)
                fps.AddFingerprint(sub_fp)
                p_sim_fps.append(sim_fp)

                p_smols.append(smol)
                p_smileses.append(smiles)

                p_links.append(links[mid])

    logging.info("Finished generating the chemical libraries")

    database = {
        "structure_smiles": p_smileses,
        "structure_wid": p_links,
        "structure_sim_fps": p_sim_fps,
        "structure_sim_h_fps": p_sim_h_fps,
        "structure_library": library.Serialize(),
        "structure_library_h": library_h.Serialize(),
        "structure_id": {i[1]: i[0] for i in enumerate(p_links)},
    }

    with open(path / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)
    logging.info("Finished dumping")

    smols_and_wids = list(zip(p_smols, p_links))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        chunks = [
            smols_and_wids[i : i + 1000] for i in range(0, len(smols_and_wids), 1000)
        ]
        sdf_blocks_list = list(executor.map(process_smol_and_wid, chunks))
        sdf_blocks = [block for sublist in sdf_blocks_list for block in sublist]
        write_mols_to_sdf(path, sdf_blocks)

    logging.info("Finished exporting")


if __name__ == "__main__":
    run(Path("data"))
