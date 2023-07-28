#!/usr/bin/env python3
import base64

import streamlit as st

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from processing_common import load_all_data

st.set_page_config(page_title="LOTUS taxon search", page_icon=":lotus:", layout="wide")


params = st.experimental_get_query_params()


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
    return load_all_data()

db = load_data()


wid = None
if "id" in params:
    if len(params["id"]) > 0:
        try:
            wid = int(params["id"][0])
        except Exception as e:
            wid = None


taxa = db["taxa"]
compounds = db["compounds"]
t2c = db["t2c"]
c2t = db["c2t"]
ccount = db["ccount"]
st.write("We have **{}** taxa 🍀 and **{}** compounds 🫧 and **{}** couples 🏩.".format(len(taxa), len(compounds), ccount))

if "taxon" not in st.session_state and wid is not None:
    if wid in db["taxa_name"]:
        st.session_state["taxon"] = db["taxa_name"][wid]

query = st.text_input("Enter the name or part of the name of your taxa of interest", key="taxon")
matches = []
query = query.lower()
if query != "" and len(query) > 3:
    for taxon in taxa.keys():
        if query in taxon.lower():
            matches.append(taxon)
else:
    st.write("Please enter a query longer than 3 characters.")
matches = sorted(matches)
for match in matches:
    for i in taxa[match]:
        if i in t2c:
            with st.expander(f"{match} - {len(t2c[i])} compound(s)", expanded=False):
                st.markdown(f"[Wikidata page of {match}](https://www.wikidata.org/entity/Q{i})")
                cs = st.columns(3)
                for idx, j in enumerate(t2c[i]):
                    with cs[idx%3]:
                        link_svg(f"/Molecule_taxa?id={j}&type=structure", molecule_svg(Chem.MolFromSmiles(compounds[j])))
                        st.markdown(f"[Wikidata page of compound](https://www.wikidata.org/entity/Q{j})")
                        taxa_count = len(c2t[j])
                        if taxa_count == 1:
                            st.markdown(f"[Compound only found in this taxon](/Molecule_taxa?id={j}&type=structure)")
                        else:
                            st.markdown(f"[Found in {taxa_count} taxa](/Molecule_taxa?id={j}&type=structure)")
                        st.divider()
