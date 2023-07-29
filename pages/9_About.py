import streamlit as st
from ui_common import on_all_pages

st.set_page_config(page_title="About LOTUS", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)
on_all_pages()


@st.cache_data(ttl=3600)
def readme():
    return "".join([line for line in open("README.md").readlines() if not line.startswith(" ")])


st.markdown(readme())
