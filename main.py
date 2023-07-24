#!/usr/bin/env python3
import base64

import streamlit as st
import csv

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def molecule_svg(mol):
    d2d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def render_svg(svg):
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)

def link_svg(link, svg):
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<a href="%s"><img src="data:image/svg+xml;base64,%s"/></a>' % (link,b64)
    st.write(html, unsafe_allow_html=True)


@st.cache_resource(ttl=3600)
def load_data():
    taxa = {}
    compounds = {}
    t2c = {}
    c2t = {}

    with open("data/couples.csv", "r") as f:
        ccount = 0
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            c, t = x
            ic = int(c)
            it = int(t)
            if it not in t2c:
                t2c[it] = set()
            if ic not in c2t:
                c2t[ic] = set()
            t2c[it].add(ic)
            c2t[ic].add(it)
            ccount += 1
    with open("data/taxa.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            t, i = x
            if i not in taxa:
                taxa[i] = set()
            taxa[i].add(int(t))


    with open("data/smiles.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)
        for x in reader:
            c, smi, cano = x
            if smi == "":
                smi = cano
            compounds[int(c)] = smi
    return taxa, compounds, t2c, c2t, ccount


taxa, compounds, t2c, c2t, ccount = load_data()

st.write("We have **{}** taxa 🍀 and **{}** compounds 🫧 and **{}** couples 🏩.".format(len(taxa), len(compounds), ccount))


query = st.text_input("Enter the name of your taxa of interest", key="taxon")
matches = []
query = query.lower()
if query != "":
    for taxon in taxa.keys():
        if query in taxon.lower():
            matches.append(taxon)
matches = sorted(matches)
for match in matches:
    for i in taxa[match]:

        with st.expander(f"{match} - {len(t2c[i])} compound(s)", expanded=False):
            st.markdown(f"[Wikidata page of {match}](https://www.wikidata.org/entity/Q{i})")
            for j in t2c[i]:
                link_svg(f"https://www.wikidata.org/entity/Q{j}", molecule_svg(Chem.MolFromSmiles(compounds[j])))
                st.write("Compound found in **{}** taxa.".format(len(c2t[j])))
                st.divider()
