#!/usr/bin/env python3
import base64

import streamlit as st

from processing_common import load_all_data
from ui_common import on_all_pages

st.set_page_config(page_title="LOTUS taxon search", page_icon=":lotus:", layout="wide")
on_all_pages()

params = st.experimental_get_query_params()

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
if query != "" and len(query) > 2:
    for taxon in taxa.keys():
        if query in taxon.lower():
            matches.append(taxon)
else:
    st.write("Please enter a query longer than 2 characters.")
matches = sorted(matches)

for match in matches:
    for i in taxa[match]:
        if i in db["t2c"]:
            matching_compounds = set(db["t2c"][i])
        else:
            matching_compounds = set()
        if i in db["taxa_childs"]:
            for parent in db["taxa_childs"][i]:
                if parent in db["t2c"]:
                    for compound in db["t2c"][parent]:
                        matching_compounds.add(compound)
        st.markdown(f"[{match} - {len(matching_compounds)} compound(s)](/taxon?wid={i}&type=taxon)")
