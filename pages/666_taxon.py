import streamlit as st
import streamlit.components.v1 as components
from chemistry_helpers import molecule_svg

from chemistry_helpers import molecule_svg
from processing_common import load_all_data
from ui_common import link_svg, on_all_pages

params = st.experimental_get_query_params()

st.set_page_config(page_title="LOTUS Structure Info", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)
on_all_pages()


@st.cache_resource(ttl=3600)
def load_data():
    return load_all_data()


db = load_data()

wid = None
try:
    if "wid" in params and "type" in params:
        if params["type"][0] == "taxon":
            if len(params["wid"]) > 0:
                wid = int(params["wid"][0])
except Exception as e:
    wid = None

if id is not None and wid is not None:
    match = db["taxa_name"][wid]
    st.header(match)
    if wid in db["t2c"]:
        matching_compounds = set(db["t2c"][wid])
    else:
        matching_compounds = set()
    if wid in db["taxa_childs"]:
        for parent in db["taxa_childs"][wid]:
            if parent in db["t2c"]:
                for compound in db["t2c"][parent]:
                    matching_compounds.add(compound)

    st.markdown(f"[Wikidata page of {match}](https://www.wikidata.org/entity/Q{wid})")
    cs = st.columns(3)
    for idx, j in enumerate(list(matching_compounds)):
        with cs[idx % 3]:
            link_svg(f"/molecule?id={j}&type=structure", molecule_svg(db["compounds"][j]))
            st.markdown(f"[Wikidata page of compound](https://www.wikidata.org/entity/Q{j})")
            taxa_count = len(db["c2t"][j])
            if taxa_count == 1:
                st.markdown(f"[Compound only found in this taxon](/molecule?id={j}&type=structure)")
            else:
                st.markdown(f"[Found in {taxa_count} taxa](/molecule?id={j}&type=structure)")
            st.divider()
