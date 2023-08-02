import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st
from sqlalchemy import text

from chemistry import molecule_svg

DATA_PATH = Path("../data")
CONTRIBUTIONS_PATH = Path("../contributions")

st.set_page_config(page_title="Lotus admin - contributions", page_icon=":lotus:", layout="wide",
                   initial_sidebar_state="auto", menu_items=None)


def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    """
    if os.path.exists(path_to_file):
        timestamp = os.stat(path_to_file).st_ctime
        return datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')
    else:
        return 'File does not exist'


conn = st.experimental_connection("contributions_db", type="sql")

# Insert some data with conn.session.
with conn.session as s:
    s.execute(text(
        'CREATE TABLE IF NOT EXISTS single_contribution (filename TEXT, molecule TEXT, organism TEXT, selected_organism TEXT, reference TEXT, email TEXT, processed bool, rejected bool, accepted bool, verified bool);'))
    # for testing
    s.execute(text('DELETE FROM single_contribution'))
    s.commit()


def load_files_into_db():
    for file in CONTRIBUTIONS_PATH.glob("*.json"):
        with open(file, "r") as f:
            content = json.load(f)
        with conn.session as s:
            query = text("SELECT COUNT(*) FROM single_contribution WHERE filename = :filename")
            count = s.execute(query, {"filename": file.name}).fetchone()[0]
            if count == 0:
                query = text(
                    "INSERT INTO single_contribution (filename, molecule, organism, selected_organism, reference, email, processed, rejected, accepted, verified) VALUES (:filename, :molecule, :organism, :selected_organism, :reference, :email, :processed, :rejected, :accepted, :verified)")
                s.execute(query,
                          {"filename": file.name, "molecule": content["molecule"], "organism": content["organism"],
                           "selected_organism": content["selected_organism"], "reference": content["reference"],
                           "email": content["email"], "processed": False, "rejected": False, "accepted": False,
                           "verified": False})
            s.commit()


load_files_into_db()
# st.header("Single contributions")
singles_df = conn.query('select * from single_contribution', ttl=0)

singles_df["structure"] = singles_df["molecule"].apply(lambda x: "data:image/svg+xml;utf8," + molecule_svg(x))
edited_df = st.data_editor(
    singles_df[["organism", "selected_organism", "molecule", "structure", "reference", "email", "processed",
                "rejected", "accepted", "verified"]],
    num_rows="dynamic",
    column_config={
        "structure": st.column_config.ImageColumn(
            "Structure", help="Structure"
        ),
        "molecule": st.column_config.TextColumn(
            "Molecule", help="Molecule", width=50
        ),
        "processed": st.column_config.CheckboxColumn(
            "Processed",
            help="Processed",
            default=False,
        ),
        "rejected": st.column_config.CheckboxColumn(
            "rejected",
            help="rejected",
            default=False,
        ),
        "accepted": st.column_config.CheckboxColumn(
            "accepted",
            help="accepted",
            default=False,
        ),
        "verified": st.column_config.CheckboxColumn(
            "verified",
            help="verified",
            default=False,
        )
    },
)
#st.dataframe(edited_df)
st.header("Multiple contributions")
for file in CONTRIBUTIONS_PATH.glob("*.xsv"):
    txt_file = file.with_suffix(".txt")

    if os.path.isfile(txt_file):
        with open(txt_file, "r") as f:
            name = f"`{f.read()}`"
    else:
        name = "Unnamed contributor"

    name += f" - {creation_date(file)}"

    with st.expander(name):

        df = pd.read_csv(file, sep='\t')
        st.dataframe(df)
