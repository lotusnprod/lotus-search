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
 
"""
)
