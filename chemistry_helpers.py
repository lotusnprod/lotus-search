import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

@st.cache_data
def mol_from_smiles(smiles):
    return Chem.MolFromSmiles(smiles)

@st.cache_data
def molecule_svg(smiles, width=250):
    mol = mol_from_smiles(smiles)
    d2d = rdMolDraw2D.MolDraw2DSVG(width, width)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()




@st.cache_data
def molecule_png(smiles, width=600):
    mol = mol_from_smiles(smiles)
    return Chem.Draw.MolToImage(mol, size=(width, width))
