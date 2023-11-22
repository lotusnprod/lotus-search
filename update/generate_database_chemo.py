#!/usr/bin/env python3
import csv
import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

from processing_common import fingerprint, standardize


def process_smiles(inp):
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
        print("Failed to process", smiles)
        return None


def run(root: Path) -> None:
    smileses = []
    links = []
    with open("./data/smiles.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            c, smi, cano = x
            if smi == "":
                smi = cano
            smileses.append(smi)
            links.append(int(c))

    max_workers = multiprocessing.cpu_count()

    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()

    library = rdSubstructLibrary.SubstructLibrary(mols, fps)

    mols_h = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps_h = rdSubstructLibrary.PatternHolder()

    library_h = rdSubstructLibrary.SubstructLibrary(mols_h, fps_h)

    p_smileses = []
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

                mols.AddSmiles(smiles_clean)
                fps.AddFingerprint(sub_fp)
                p_sim_fps.append(sim_fp)

                p_smileses.append(smiles)

                p_links.append(links[mid])

    t2c = {}
    c2t = {}
    tc2r = {}

    with open("./data/couples.csv", "r") as f:
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
            tc2r[(it, ic)] = r
    print("Finished generating")

    database = {
        "structure_smiles": p_smileses,
        "structure_wid": p_links,
        "structure_sim_fps": p_sim_fps,
        "structure_sim_h_fps": p_sim_h_fps,
        "structure_library": library.Serialize(),
        "structure_library_h": library_h.Serialize(),
        "t2c": t2c,
        "c2t": c2t,
        "tc2r": tc2r,
    }

    with open(root / "database_chemo.pkl", "wb") as f:
        pickle.dump(database, f)
    print("Finished dumping")


if __name__ == "__main__":
    run(Path("data"))
