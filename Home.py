import streamlit as st

from ui_common import on_all_pages

st.set_page_config(
    page_title="LOTUS",
    page_icon="👋",
)

on_all_pages()

st.write("# Welcome to LOTUS! 👋")

st.sidebar.success("Select what you want to search for above")

st.markdown(
    """
    We explore new ways to share knowledge in Natural Products research. 
    
    This is LOTUS, it aims at sharing an online resource for NP 
    occurrences. The LOTUS data is available under [CC0 license](https://creativecommons.org/publicdomain/zero/1.0/). 
    It is hosted in parallel at https://www.wikidata.org/ and 
    https://search.nprod.net/.
    
    The Wikidata version allows for community curation and addition
     of novel data and this search engine is updated daily with new
     data. 
     
     You can currently:
     
     [Search for a taxon](/Taxon_search)
     
     and
     
     [Search for a structure](/Structure_search)
     
     (or use the buttons on the left, click on the arrow on the top left of your screen if you are on a phone)
"""
)
