import json
import os
import random
import time
import uuid
from io import StringIO
from pathlib import Path

import streamlit as st
from bjonnh_streamlit_ketcher import st_ketcher

from chemistry_helpers import molecule_svg
from ui_common import data_model, on_all_pages

st.set_page_config(page_title="LOTUS Contribution", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)
on_all_pages()

CONTRIBUTIONS_PATH = Path("./contributions/")


dm = data_model()


def captcha_solver():
    if "captcha-code" not in st.session_state:
        st.session_state["captcha-code"] = ''.join(random.choices("CONPS", k=5))
    image = molecule_svg(st.session_state['captcha-code'])
    col1, col2 = st.columns(2)
    col1.image(image)
    captcha_user = col2.text_input('Enter the molecule SMILES string (without the hydrogens!) and click on Verify')
    if st.button("Verify"):
        if captcha_user != st.session_state["captcha-code"] and captcha_user != st.session_state["captcha-code"][::-1]:
            st.session_state["captcha"] = False
            st.error("Sorry that's invalid")
        else:
            st.session_state["captcha"] = True
            st.experimental_rerun()


def ui_email() -> str | None:
    st.subheader("Your email (optional)")
    st.write(
        "We will store it for the purpose of contacting you if we have questions it is not shared to third parties")
    return st.text_input("Email")


st.header("Awesome!")
st.write("""That means you will be able to help us improve the data quality of LOTUS. We trust that you have checked that
        it is not already present in the database and that your molecule is drawn correctly.""")
st.write("""It is preferable for you to add the molecules directly on Wikidata, we are working on a new
tutorial to make that easy for you. In the meantime, please fill the form below and we will add it for you!""")

disk = sum(os.path.getsize(CONTRIBUTIONS_PATH/f)+1 for f in os.listdir(CONTRIBUTIONS_PATH) if os.path.isfile(CONTRIBUTIONS_PATH/f))
nf = len(os.listdir(CONTRIBUTIONS_PATH))
if disk > 100_000_000 or nf > 1000:
    st.write("Unfortunately we already have too many contributions, please try again later")
    st.write(f"Error code: S{disk}F{nf}")
    st.stop()

if "captcha" not in st.session_state or st.session_state["captcha"] is False:
    st.error("Please solve the captcha below to prove you are a human")
    captcha_solver()
else:
    email = ui_email()

    select = st.selectbox("Do you want to add a single compound or multiple?", ["Single", "Multiple"])
    if select == "Single":
        st.header("1. Draw your compound")
        molecule = st_ketcher()
        if molecule:
            st.header("2. Give the name of the organism")
            organism = st.text_input(
                "Organism (please use verify and select from the list, pick the most specific name you can find)")
            selected_organism = None
            if st.button("Verify") or "selected_organism" in st.session_state:
                results = dm.get_taxa_with_name_containing(organism)
                names = sorted([dm.get_taxon_name_from_wid(match) for match in results])
                selected_organism = st.selectbox("Select the most specific name", names, key="selected_organism")

            if organism:
                st.header("3. Give the reference to the paper (please ideally a DOI)")
                reference = st.text_input("Reference (DOI ideally)")

                if reference:
                    if len(molecule) > 8192:
                        st.error("Your molecule is too big, we do not handle molecules this size yet.")
                    elif len(organism) > 100:
                        st.error("Your organism name is too long, please use a shorter name.")
                    elif len(reference) > 512:
                        st.error("Your reference is way too long, please use a doi or a short reference.")
                    elif len(email) > 64 or ("@" not in email and email != ""):
                        st.error("Your email is invalid, please use a valid email address or none at all")
                    if st.button("Submit"):
                        js = json.dumps({
                            "molecule": molecule,
                            "organism": organism,
                            "selected_organism": selected_organism,
                            "reference": reference,
                            "email": email
                        })
                        with open(f"contributions/{time.time()}_{uuid.uuid4()}.json", "w") as f:
                            f.write(js)
                        st.success("Thank you for your contribution, we will review and eventually contact you!")

                        st.write("If you are curious, this is the only information we store about your request:")
                        st.json(js)

    else:
        st.write("You are even more awesome!")
        st.write(
            "So if you have way too many compounds to add, you can send us directly a TSV or CSV file with the following columns:")
        st.markdown("**smiles** , **organism** , **reference**")
        st.write(
            "Please make sure that it is formated properly and that the values are quoted (even if they don't contain commas)")

        uploaded_file = st.file_uploader("Upload your file", type=["tsv", "csv"])
        if uploaded_file is not None:

            if st.button("Send"):
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                basename = f"contributions/{time.time()}_{uuid.uuid4()}"
                with open(f"{basename}.xsv", "w") as f:
                    f.write(stringio.read())
                if email != "" and email is not None and "@" in email and len(email)< 64:
                    with open(f"{basename}.txt", "w") as f:
                        f.write(email)
                st.success("Thanks a lot for your contributions, we will eventually contact you if we have questions")
