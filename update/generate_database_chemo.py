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

from chemistry_helpers import process_smiles, write_mols_to_sdf

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run(path: Path) -> None:
    smileses = []
    links = []
    with open(path / "structures.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        # for x in islice(reader, 1000):
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

    sdf_blocks = []
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
                    nid,
                    smiles,
                    smol,
                    smiles_clean,
                    mol_block,
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
                sdf_blocks.append((links[nid], mol_block))
                fps.AddFingerprint(sub_fp)
                p_sim_fps.append(sim_fp)

                p_smols.append(smol)
                p_smileses.append(smiles)

                p_links.append(links[nid])

    logging.info("Finished generating the chemical libraries")

    with open(path / "smiles_processed.csv", "w") as f:
        # Write a csv with header, structure_id and structure_smiles
        # from the two arrays p_links and p_smileses  respectively
        csv_file = csv.writer(f)
        csv_file.writerow(["structure", "structure_smiles"])
        csv_file.writerows(zip(p_links, p_smileses))

    database = {
        "structure_wid": p_links,
        # TODO add blocks if needed
        "structure_sim_fps": p_sim_fps,
        "structure_sim_h_fps": p_sim_h_fps,
        "structure_library": library.Serialize(),
        "structure_library_h": library_h.Serialize(),
        "structure_id": {i[1]: i[0] for i in enumerate(p_links)},
    }

    with open(path / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)
    logging.info("Finished dumping")

    logging.info("Exporting SDF")
    write_mols_to_sdf(path, sdf_blocks)

    logging.info("Finished exporting")


if __name__ == "__main__":
    run(Path("data"))
