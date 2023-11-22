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
    new_lib_h = rdSubstructLibrary.SubstructLibrary()
    with io.BytesIO(data["structure_library"]) as i:
        new_lib.InitFromStream(i)
    with io.BytesIO(data["structure_library_h"]) as i:
        new_lib_h.InitFromStream(i)
    data["structure_library"] = new_lib
    data["structure_library_h"] = new_lib_h

    data["structure_id"] = {i[1]: i[0] for i in enumerate(data["structure_wid"])}
    return data
