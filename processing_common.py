import csv
import pickle

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, rdSubstructLibrary
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
uncharger = rdMolStandardize.Uncharger()


def standardize(mol):
    clean_mol = rdMolStandardize.Cleanup(mol)
    bigger_clean = rdMolStandardize.FragmentParent(clean_mol)
    bigger_clean = uncharger.uncharge(bigger_clean)
    return bigger_clean


def fingerprint(mol):
    return fpgen.GetFingerprint(mol)


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

def load_all_data():
    taxa = {}
    taxa_childs = {}
    compounds = {}
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
    with open("./data/taxa.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            taxon_id, taxon_name, parent_id, parent_name = x
            if taxon_name not in taxa:
                taxa[taxon_name] = set()
            if parent_name not in taxa:
                taxa[parent_name] = set()
            if parent_id not in taxa_childs:
                taxa_childs[parent_id] = set()
            taxa[taxon_name].add(int(taxon_id))
            taxa[parent_name].add(int(parent_id))
            taxa_childs[parent_id].add(int(taxon_id))

    with open("./data/smiles.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            c, smi, cano = x
            if smi == "":
                smi = cano
            compounds[int(c)] = smi

    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()
    with open("./data/structure_fps.pkl", "rb") as f:
        data = pickle.load(f)
    for idx in range(len(data["sub_fps"])):
        mols.AddSmiles(data["smileses_clean"][idx])
        fps.AddFingerprint(data["sub_fps"][idx])

    data["library"] = rdSubstructLibrary.SubstructLibrary(mols, fps)
    data["taxa"] = taxa
    data["taxa_name"] = {}
    for taxa_name in data["taxa"].keys():
        for i in data["taxa"][taxa_name]:
            data["taxa_name"][i] = taxa_name
    data["compounds"] = compounds
    data["t2c"] = t2c
    data["c2t"] = c2t
    data["tc2r"] = tc2r
    data["taxa_childs"] = taxa_childs
    data["ccount"] = ccount
    return data
