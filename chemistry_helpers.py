from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


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
