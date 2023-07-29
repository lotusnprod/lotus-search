import streamlit as st
import streamlit.components.v1 as components

from chemistry_helpers import molecule_svg
from processing_common import load_all_data

params = st.experimental_get_query_params()

st.set_page_config(page_title="LOTUS Structure Info", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)

@st.cache_resource(ttl=3600)
def load_data():
    return load_all_data()

db = load_data()

id = None
wid = None
try:
    if "id" in params and "type" in params:
        if params["type"][0] == "structure":
            if len(params["id"]) > 0:
                wid = int(params["id"][0])
                id = db["links"].index(str(wid))
except Exception as e:
    id = None
    wid = None

if id is not None and wid is not None:
    m = db["smileses"][id]
    st.image(molecule_svg(m, width=250), use_column_width="always", width=250)
    st.header("Species")
    name_id_list = []
    for t in db["c2t"][wid]:
        name = db["taxa_name"][t]
        name_id_list.append([name, t])
    name_id_list = sorted(name_id_list, key=lambda x: x[0])
    for name, wid_tax in name_id_list:
        st.markdown(f"[{name}](/Taxon_search?id={wid_tax})")
    st.header("Wikidata")
    st.write("You can find a lot more information on the page of this compound on Wikidata.")
    st.markdown(f"[Go to the wikidata page](http://www.wikidata.org/entity/Q{wid})")
    components.iframe(f"https://www.wikidata.org/entity/Q{wid}", scrolling=True, height=500)
    if id in db["c2t"]:
        taxa_count = len(db["c2t"][id])
        st.markdown(f"Found in {taxa_count} taxa")