from rdkit import Chem
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


def standardize(mol):
    clean_mol = rdMolStandardize.Cleanup(mol)
    bigger_clean = rdMolStandardize.FragmentParent(clean_mol)
    bigger_clean = uncharger.uncharge(bigger_clean)
    return bigger_clean


fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
uncharger = rdMolStandardize.Uncharger()
