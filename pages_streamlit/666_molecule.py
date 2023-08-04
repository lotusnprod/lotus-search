import streamlit as st
import streamlit.components.v1 as components

from chemistry_helpers import molecule_svg
from config import TTL_DATA_CACHE
from model import DataModel
from processing_common import load_all_data
from ui_common import get_url_parameter, on_all_pages

params = st.experimental_get_query_params()

st.set_page_config(page_title="LOTUS Structure Info", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)
on_all_pages()


@st.cache_resource(ttl=TTL_DATA_CACHE)
def data_model():
    return DataModel(load_all_data())


dm = data_model()

wid = get_url_parameter("id", "structure")

if wid is not None:
    st.header(f"Q{wid}")

    m = dm.get_compound_smiles_from_wid(wid)  ## TODO Switch to mol
    st.image(molecule_svg(m, width=250), use_column_width="always", width=250)
    c1, c2 = st.columns(2)
    c1.markdown(f"[Load in structure editor](/Structure_search?id={wid}&type=structure)")
    c2.markdown(f"[Go to the Wikidata page](http://www.wikidata.org/entity/Q{wid})")
    st.header("Found in species")
    name_id_list = []
    for t in dm.get_taxa_containing_compound(wid):
        name = dm.get_taxon_name_from_wid(t)
        name_id_list.append([name, t])
    name_id_list = sorted(name_id_list, key=lambda x: x[0])
    for name, wid_tax in name_id_list:
        st.markdown(f"[{name}](/taxon?wid={wid_tax}&type=taxon)")
    st.header("Wikidata")
    st.write("You can find a lot more information on the page of this compound on Wikidata.")
    st.markdown(f"[Go to the Wikidata page](http://www.wikidata.org/entity/Q{wid})")
    components.iframe(f"https://www.wikidata.org/entity/Q{wid}", scrolling=True, height=500)
