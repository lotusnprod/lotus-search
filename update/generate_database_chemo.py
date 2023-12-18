#!/usr/bin/env python3
import csv
import logging
import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdSubstructLibrary

from chemistry_helpers import fingerprint, standardize

RDLogger.DisableLog("rdApp.*")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def process_smol_and_wid(smol_and_wid):
    sdf_blocks = []
    for smol, wid in smol_and_wid:
        sdf_blocks.append((wid, Chem.MolToMolBlock(smol)))
    return sdf_blocks


def process_smiles(inp):
    smiles = "Input to process_smiles is invalid"
    try:
        nid, smiles = inp
        mol = Chem.MolFromSmiles(smiles)
        smol = standardize(mol)
        smiles_clean = Chem.MolToSmiles(smol)
        if smol is not None:
            sim_fp = fingerprint(smol)
            sub_fp = Chem.PatternFingerprint(smol)
            smol_h = Chem.AddHs(smol)
            sim_fp_h = fingerprint(smol_h)
            sub_fp_h = Chem.PatternFingerprint(smol_h)
            return (
                nid,
                smiles,
                smol,
                smiles_clean,
                sim_fp,
                sub_fp,
                smol_h.ToBinary(),
                sim_fp_h,
                sub_fp_h,
            )
        else:
            return None
    except:
        logging.error(f"Failed to process: {smiles}")
        return None


def write_mols_to_sdf(path: Path, sdf_blocks):
    with Chem.SDWriter(str(path / "lotus.sdf")) as w:
        for wid, sdf_block in sdf_blocks:
            mol = Chem.MolFromMolBlock(sdf_block)
            if mol:
                mol.SetProp("WID", str(wid))
                w.write(mol)


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

                mols_h.AddMol(Chem.Mol(mol_h))
                fps_h.AddFingerprint(sub_fp_h)
                p_sim_h_fps.append(sim_fp_h)

                smis.AddSmiles(smiles_clean)
                fps.AddFingerprint(sub_fp)
                p_sim_fps.append(sim_fp)

                p_smols.append(smol)
                p_smileses.append(smiles)

                p_links.append(links[mid])

    t2c: dict[int, set[int]] = {}
    c2t: dict[int, set[int]] = {}
    tc2r: dict[tuple, set[str]] = {}

    logging.info("Finished generating the chemical libraries")

    with open(path / "couples.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            c, t, r = x
            ic = int(c)
            it = int(t)
            if it not in t2c:
                t2c[it] = set()
            if ic not in c2t:
                c2t[ic] = set()
            if (it, ic) not in tc2r:
                tc2r[(it, ic)] = set()
            t2c[it].add(ic)
            c2t[ic].add(it)
            tc2r[(it, ic)] = set(r)
    logging.info("Finished generating couples")

    database = {
        "structure_smiles": p_smileses,
        "structure_wid": p_links,
        "structure_sim_fps": p_sim_fps,
        "structure_sim_h_fps": p_sim_h_fps,
        "structure_library": library.Serialize(),
        "structure_library_h": library_h.Serialize(),
        "structure_id": {i[1]: i[0] for i in enumerate(p_links)},
        "t2c": t2c,
        "c2t": c2t,
        "tc2r": tc2r,
    }

    smols_and_wids = list(zip(p_smols, p_links))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        chunks = [
            smols_and_wids[i : i + 1000] for i in range(0, len(smols_and_wids), 1000)
        ]
        sdf_blocks_list = list(executor.map(process_smol_and_wid, chunks))
        sdf_blocks = [block for sublist in sdf_blocks_list for block in sublist]
        write_mols_to_sdf(path, sdf_blocks)

    logging.info("Finished exporting")

    with open(path / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)
    logging.info("Finished dumping")


if __name__ == "__main__":
    run(Path("data"))
