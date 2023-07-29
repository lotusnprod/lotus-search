import csv
import pickle

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, rdSubstructLibrary
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

# Memory is cheap!
def load_all_data():
    compounds = {}

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

    with open("./data/database.pkl", "rb") as f:
        data = pickle.load(f)
    for idx in range(len(data["sub_fps"])):
        mols.AddSmiles(data["smileses_clean"][idx])
        fps.AddFingerprint(data["sub_fps"][idx])

    data["library"] = rdSubstructLibrary.SubstructLibrary(mols, fps)

    data["taxa_name"] = {}
    for taxa_name in data["taxa"].keys():
        for i in data["taxa"][taxa_name]:
            data["taxa_name"][i] = taxa_name
    data["compounds"] = compounds
    return data
