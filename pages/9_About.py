import streamlit as st

@st.cache_data(ttl=3600)
def readme():
    return "".join([line for line in open("README.md").readlines() if not line.startswith(" ")])


st.markdown(readme())
