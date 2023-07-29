#!/usr/bin/env python3
import csv
import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem

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
            return (nid, smiles, smiles_clean, sim_fp, sub_fp)
        else:
            return None
    except:
        print("Failed to process", smiles)
        return None


def run(root: Path) -> None:
    with open(root / "smiles.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        stuff = [x for x in reader]
        links, isomeric_smileses, canonical_smileses = zip(*stuff)
        smileses = []
        for i in range(len(isomeric_smileses)):
            if isomeric_smileses[i] == "":
                smileses.append(canonical_smileses[i])
            else:
                smileses.append(isomeric_smileses[i])

    max_workers = multiprocessing.cpu_count()

    p_smileses = []
    p_smileses_clean = []
    p_sim_fps = []
    p_sub_fps = []
    p_links = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(process_smiles, enumerate(smileses), chunksize=1000)
        for result in results:
            if result is not None:
                mid, smiles, smiles_clean, sim_fp, sub_fp = result
                p_smileses.append(smiles)
                p_smileses_clean.append(smiles_clean)
                p_sim_fps.append(sim_fp)
                p_sub_fps.append(sub_fp)
                p_links.append(links[mid].replace("http://www.wikidata.org/entity/", ""))

    taxa = {}
    taxa_childs = {}
    taxa_parents = {}
    with open(root / "taxa.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            taxon_id, taxon_name, parent_id, parent_name = x
            taxon_id = int(taxon_id)
            parent_id = int(parent_id)
            if taxon_name not in taxa:
                taxa[taxon_name] = set()
            if parent_name not in taxa:
                taxa[parent_name] = set()
            if parent_id not in taxa_childs:
                taxa_childs[parent_id] = set()
            if taxon_id not in taxa_parents:
                taxa_parents[taxon_id] = set()
            taxa[taxon_name].add(taxon_id)
            taxa[parent_name].add(parent_id)
            taxa_childs[parent_id].add(taxon_id)
            taxa_parents[taxon_id].add(parent_id)

    t2c = {}
    c2t = {}
    tc2r = {}

    with open("./data/couples.csv", "r") as f:
        ccount = 0
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
                ccount += 1
            t2c[it].add(ic)
            c2t[ic].add(it)
            tc2r[(it, ic)] = r

    database = {
        "smileses": p_smileses,
        "smileses_clean": p_smileses_clean,
        "links": p_links,
        "sim_fps": p_sim_fps,
        "sub_fps": p_sub_fps,
        "taxa": taxa,
        "taxa_childs": taxa_childs,
        "taxa_parents": taxa_parents,
        "t2c": t2c,
        "c2t": c2t,
        "tc2r": tc2r,
        "ccount": ccount
    }

    pickle.dump(database, open(root / "database.pkl", "wb"))
