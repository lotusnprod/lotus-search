#!/usr/bin/env python3
import base64

import streamlit as st

from config import TTL_DATA_CACHE
from model import DataModel
from processing_common import load_all_data
from ui_common import data_model, get_url_parameter, on_all_pages

st.set_page_config(page_title="LOTUS taxon search", page_icon=":lotus:", layout="wide")
on_all_pages()

dm = data_model()

wid = get_url_parameter("id", "taxon")

st.write("We have **{}** taxa 🍀 and **{}** compounds 🫧 and **{}** couples 🏩.".format(dm.num_taxa(),
                                                                                     dm.num_compounds(),
                                                                                     dm.num_couples()))

if "taxon" not in st.session_state and wid is not None:
    name = dm.get_taxon_name_from_wid(wid)
    if name is not None:
        st.session_state["taxon"] = name

query = st.text_input("Enter the name or part of the name of your taxa of interest", key="taxon")
matches = []

if query != "" and len(query) > 2:
    matches.extend(dm.get_taxa_with_name_containing(query))
else:
    st.write("Please enter a query longer than 2 characters.")
matches = sorted(matches)

for match in matches:
    name = dm.get_taxon_name_from_wid(match)
    matching_compounds = dm.get_compounds_of_taxon(match)

    st.markdown(f"[{name} - {len(matching_compounds)} compound(s)](/taxon?wid={match}&type=taxon)")
