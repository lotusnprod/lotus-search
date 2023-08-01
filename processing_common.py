import io
import pickle

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
    with open("./data/database.pkl", "rb") as f:
        data = pickle.load(f)
    new_lib = rdSubstructLibrary.SubstructLibrary()
    with io.BytesIO(data["compound_library"]) as i:
        new_lib.InitFromStream(i)
    data["compound_library"] = new_lib

    data["compound_id"] = {i[1]: i[0] for i in enumerate(data["compound_wid"])}
    return data
