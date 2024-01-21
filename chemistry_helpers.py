import logging
from mordred import Calculator, descriptors
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize


def molecule_svg(smiles: str, molecule: str | None, width: int = 250):
    mol = Chem.MolFromSmiles(smiles)
    d2d = rdMolDraw2D.MolDraw2DSVG(width, width)
    if molecule is not None:
        explicit_h = molecule is None or "[H]" in molecule
        p = Chem.SmilesParserParams()
        p.removeHs = not explicit_h
        if explicit_h:
            mol = Chem.AddHs(mol)
        highlight_atoms = mol.GetSubstructMatch(Chem.MolFromSmiles(molecule, p))
        draw_options = d2d.drawOptions()
        draw_options.setHighlightColour((0.1, 0.9, 0.9, 0.8))
        d2d.DrawMolecule(mol, highlightAtoms=highlight_atoms)
    else:
        d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def fingerprint(mol):
    return fpgen.GetFingerprint(mol)


# from https://github.com/mordred-descriptor/mordred/tree/develop/examples
def get_mol_descriptors_mordred(mol):
    return calc(mol).drop_missing().asdict()


# from https://greglandrum.github.io/rdkit-blog/posts/2022-12-23-descriptor-tutorial.html
def get_mol_descriptors_rdkit(mol):
    """ calculate the full list of descriptors for a molecule
    
        missingVal is used if the descriptor cannot be calculated
    """
    res = {}
    for nm,fn in Descriptors._descList:
        # some of the descriptor fucntions can throw errors if they fail, catch those here:
        try:
            val = fn(mol)
        except:
            # print the error message:
            import traceback
            traceback.print_exc()
            # and set the descriptor value to None
            val = None
        res[nm] = val
    return res


def standardize(mol):
    clean_mol = rdMolStandardize.Cleanup(mol)
    bigger_clean = rdMolStandardize.FragmentParent(clean_mol)
    bigger_clean = uncharger.uncharge(bigger_clean)
    return bigger_clean


def process_smiles(inp):
    smiles = "Input to process_smiles is invalid"
    try:
        nid, smiles = inp
        mol = Chem.MolFromSmiles(smiles)
        smol = standardize(mol)
        if smol is not None:
            smiles_clean = Chem.MolToSmiles(smol)
            inchi_clean = Chem.inchi.MolToInchi(smol)
            inchikey_clean = Chem.inchi.MolToInchiKey(smol)
            mol_block = Chem.MolToMolBlock(smol)
            sim_fp = fingerprint(smol)
            sub_fp = Chem.PatternFingerprint(smol)
            desc_mordred = get_mol_descriptors_mordred(smol)
            desc_rdkit = get_mol_descriptors_rdkit(smol)
            smol_h = Chem.AddHs(smol)
            sim_fp_h = fingerprint(smol_h)
            sub_fp_h = Chem.PatternFingerprint(smol_h)
            Chem.RemoveStereochemistry(smol)
            smiles_no_stereo = Chem.MolToSmiles(smol)
            inchi_no_stereo = Chem.inchi.MolToInchi(smol)
            return (
                nid,
                smiles,
                smol,
                smiles_clean,
                inchi_clean,
                inchikey_clean,
                mol_block,
                sim_fp,
                sub_fp,
                desc_mordred,
                desc_rdkit,
                smol_h.ToBinary(),
                sim_fp_h,
                sub_fp_h,
                smiles_no_stereo,
                inchi_no_stereo,
            )
        else:
            return None
    except:
        logging.error(f"Failed to process: {smiles}")
        return None


calc = Calculator(descriptors)
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
uncharger = rdMolStandardize.Uncharger()
