import streamlit as st

from chemistry_helpers import molecule_svg
from ui_common import data_model, get_url_parameter, link_svg, on_all_pages
from streamlit.runtime.scriptrunner.script_run_context import get_script_run_ctx
params = st.experimental_get_query_params()

dm = data_model()

wid = get_url_parameter("wid", "taxon")
if wid is not None:
    name = dm.get_taxon_name_from_wid(wid)

get_script_run_ctx()._set_page_config_allowed = True  ## We can we only got the parameters
st.set_page_config(page_title=f"LOTUS Taxon Info - {name}", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)
on_all_pages()


def tsv(compounds: list[int]) -> str:
    smileses = dm.get_compound_smiles_from_list_of_wid(compounds)
    return "smiles\n"+"\n".join(smileses)


if wid is not None:
    ranks = dm.get_ranks_string(wid)
    st.header(f"{name}{ranks}")

    matching_compounds = dm.get_compounds_of_taxon(wid)
    if wid in dm.db["taxonomy_direct_parents"]:
        markdown = ""
        parent_taxa = dm.db["taxonomy_direct_parents"][wid]
        tree = []
        for parent in parent_taxa:
            tree.append([parent, 1])
            if parent in dm.db["taxonomy_parents_with_distance"]:
                for relative in dm.db["taxonomy_parents_with_distance"][parent]:
                    distance = dm.db["taxonomy_parents_with_distance"][parent][relative]
                    tree.append([relative, distance])
        tree = sorted(tree, key=lambda x: x[1])
        for parent in tree:
            ranks = dm.get_ranks_string(parent[0])

            tax_name = dm.get_taxon_name_from_wid(parent[0])
            if tax_name is not None:
                markdown += f"[{tax_name}{ranks}](/taxon?wid={parent[0]}&type=taxon) > "
        st.markdown(markdown.strip("> "))
    st.divider()
    st.markdown(f"[Wikidata page of {name}](https://www.wikidata.org/entity/Q{wid})")

    matching_compounds = list(matching_compounds)
    if len(matching_compounds) > 100:
        st.warning(f"Found {len(matching_compounds)} compounds, but we are only showing you 100")
    else:
        st.write(f"Found {len(matching_compounds)} compounds")

    c1, c2 = st.columns(2)
    displayed = matching_compounds[0:100]
    data_dl = displayed
    if len(matching_compounds) >= 100:
        if c2.checkbox("I want to download them all and I understand it can be really slow"):
            data_dl = matching_compounds

    c1.download_button(f"Download {len(data_dl)} smiles as TSV", tsv(data_dl),
                       f"Compounds from {name.replace('.', '')}.tsv",
                       "text/tab-separated-values")

    cs = st.columns(3)
    for idx, j in enumerate(displayed):
        with cs[idx % 3]:
            link_svg(f"/molecule?id={j}&type=structure", molecule_svg(dm.get_compound_smiles_from_wid(j)))  ## TODO switch to mol
            st.markdown(f"[Wikidata page of compound](https://www.wikidata.org/entity/Q{j})")
            st.markdown(f"[Load in structure editor](/Structure_search?id={j}&type=structure)")
            taxa_count = dm.get_number_of_taxa_containing_compound(j)

            st.markdown(f"[Found in {taxa_count} taxa](/molecule?id={j}&type=structure)")
            st.divider()
