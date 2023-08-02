from urllib.parse import quote

import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import xml.etree.ElementTree as ET

@st.cache_data
def mol_from_smiles(smiles):
    return Chem.MolFromSmiles(smiles)


@st.cache_data
def molecule_svg(smiles, width=250):
    mol = mol_from_smiles(smiles)
    if mol is not None:
        d2d = rdMolDraw2D.MolDraw2DSVG(width, width)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        svg = d2d.GetDrawingText()
        return quote(svg)
    return ""
