#!/usr/bin/env python3
import time
import traceback

import streamlit as st
from bjonnh_streamlit_ketcher import st_ketcher
from rdkit import Chem

from chemistry_helpers import molecule_png
from processing_common import fingerprint, standardize
from ui_common import data_model, get_url_parameter, on_all_pages

start = time.time()
st.set_page_config(page_title="LOTUS Structure Search", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)
on_all_pages()


if "ignore_id" in st.session_state and st.session_state["ignore_id"]:
    wid = None
else:
    wid = get_url_parameter("id", "structure")

dm = data_model()


@st.cache_data(ttl=3600)
def search(fp: bytes) -> list[tuple[int, float]]:
    return dm.compound_search(fp)


@st.cache_data(ttl=3600)
def ss_search(fp: bytes, mol, chirality: bool) -> list[tuple[int, float]]:
    return dm.compound_search_substructure(fp, mol, chirality)


def tsv(scores):
    return dm.compound_get_tsv_from_scores(scores)


st.title("LOTUS structure search")

st.write("Choose an example, or draw a molecule below")
c1, c2, c3, c4 = st.columns(4)

amarogentin = "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O"

if wid is not None:
    st.session_state["input_query"] = dm.get_compound_smiles_from_wid(int(wid))

if "input_query" not in st.session_state:
    st.session_state["input_query"] = ""

if c1.button("Amarogentin"):
    st.experimental_set_query_params(id=None)
    st.session_state["input_query"] = amarogentin
    st.experimental_rerun()
if c2.button("Quassin"):
    st.experimental_set_query_params(id=None)
    st.session_state[
        "input_query"] = "C[C@@H]1C=C(C(=O)[C@]2([C@H]1C[C@@H]3[C@@]4([C@@H]2C(=O)C(=C([C@@H]4CC(=O)O3)C)OC)C)C)OC"
    st.experimental_rerun()
if c3.button("Absinthin"):
    st.experimental_set_query_params(id=None)
    st.session_state[
        "input_query"] = "C[C@H]1[C@@H]2CC[C@]([C@@H]3[C@H]4[C@H]5C=C([C@@]6([C@H]4C(=C3[C@H]2OC1=O)C)[C@@H]5[C@@](CC[C@@H]7[C@@H]6OC(=O)[C@H]7C)(C)O)C)(C)O"
    st.experimental_rerun()
if c4.button("Quinine"):
    st.experimental_set_query_params(id=None)
    st.session_state["input_query"] = "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O"
    st.experimental_rerun()

query = st_ketcher(st.session_state["input_query"])
c1, c2 = st.columns(2)
ss_mode = c1.checkbox("Sub-structure search")
if ss_mode:
    chirality = c2.checkbox("Respect chirality (careful in some cases it fails)")

if query:
    try:
        # Sometimes ketcher gives really invalid smiles like with theobromine
        mol = standardize(Chem.MolFromSmiles(query))

        fp = fingerprint(mol)

        if ss_mode:
            scores = ss_search(fp, mol, chirality)
        else:
            scores = search(fp)

            level = st.slider("Level (higher means tighter search)", min_value=0.0, max_value=1.0, value=0.8, step=0.01)
            scores = [score for score in scores if score[1] >= level]
        st.header("Results")

        start = time.time()

        count = len(scores)
        scores_sorted = sorted(scores, reverse=True, key=lambda x: x[1])
        st.write(f"Search time: {time.time() - start:.2f}s")
        scores = scores_sorted[0:100]

        if count >= 100:
            st.warning(
                f"There are {len(scores_sorted)} results, but we only show the top 100. Please refine your search.")

        c1, c2 = st.columns(2)

        scores_dl = scores
        if count >= 100:
            if c2.checkbox("I want to download them all and I understand it can be really slow"):
                scores_dl = scores_sorted

        c1.download_button(f"Download {len(scores_dl)} compounds as TSV", tsv(scores_dl), "results.tsv",
                           "text/tab-separated-values")

        cs = st.columns(4)
        for idx, result in enumerate(scores):
            with cs[idx % 4]:
                wid, score = result
                wid = int(wid)
                st.divider()
                m = dm.get_compound_smiles_from_wid(wid)  ## TODO Switch to mol
                st.image(molecule_png(m), use_column_width="always", width=250)

                cc1, cc2 = st.columns(2)
                cc1.markdown(f"[Wikidata page](http://www.wikidata.org/entity/Q{wid})")

                taxa_count = dm.get_number_of_taxa_containing_compound(wid)
                if taxa_count > 0:
                    cc1.markdown(f"[Found in {taxa_count} taxa](/molecule?id={wid}&type=structure)")

                if cc2.button("Load in editor", key=f"load_{wid}_{idx}"):
                    st.session_state["input_query"] = m
                    st.session_state["ignore_id"] = True
                    st.experimental_rerun()

                st.progress(score, text=f"Tanimoto similarity: {score:.2f}")

    except Exception as e:
        print(f"Got an exception with entry {query} and exception {e}")
        print(traceback.format_exc())
        st.error("Your molecule is likely invalid.")
