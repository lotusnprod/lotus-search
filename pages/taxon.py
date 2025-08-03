import math
from typing import Any

import dash_bootstrap_components as dbc
from dash_common import generate_structures_cards
from dash_config import PAGE_SIZE

import dash
from dash import Input, Output, callback, dcc
from model.model import DataModel

dm = DataModel()


def title(wid=None):
    if wid is not None:
        try:
            return f"LOTUS - {dm.get_taxon_object_from_tid(int(wid))}"
        except:
            pass
    return "LOTUS"


def tsv(structures: list[int]) -> str:
    smileses = dm.get_structure_smiles_from_list_of_sids(structures)
    return "smiles\n" + "\n".join(smileses)


dash.register_page(
    __name__,
    name="Taxon information",
    top_nav=True,
    order=-1,
    path_template="/taxon/<wid>",
    title=title,
)


@callback(
    Output("cards", "children"),
    [Input("pagination", "active_page"), Input("matching-ids", "data")],
)
def structure_cards(active_page: int, data: dict[str, Any]) -> list[dbc.Card]:
    return generate_structures_cards(active_page, data)


@callback(
    Output("download", "data"),
    [Input("btn-download", "n_clicks"), Input("matching-ids", "data")],
    prevent_initial_call=True,
)
def func(n_clicks, data):
    filename = f"structures_of_{data['taxon_name'].replace('.', '')}.tsv"
    return dict(content=tsv(data["matching_structures"]), filename=filename)


def layout(wid: int):
    if wid is None:
        return dbc.Container([])
    try:
        wid = int(wid)
    except ValueError:
        return dbc.Container([])

    taxon_name = dm.get_taxon_object_from_tid(wid)
    parent_ranks = dm.get_ranks_string(wid)

    taxonomic_info = ""

    tree = dm.get_taxonomic_tree(wid)

    if len(tree) > 0:
        markdown = ""
        for parent in tree:
            ranks = dm.get_ranks_string(parent[0])

            tax_name = dm.get_taxon_object_from_tid(parent[0])
            if tax_name is not None:
                markdown += f"[{tax_name}{ranks}](/taxon/{parent[0]}) > "
        taxonomic_info = markdown.strip("> ")

    matching_structures = dm.get_structures_of_taxon(wid)
    matching_structures.sort()
    nb_matches = len(matching_structures)

    warning = f"Found {nb_matches} structures"

    return dbc.Container([
        dcc.Store(
            id="matching-ids",
            data={
                "matching_structures": matching_structures,
                "taxon_name": taxon_name,
            },
        ),
        dbc.Row([
            dash.html.H1(f"{taxon_name}{parent_ranks}"),
            dcc.Markdown(taxonomic_info),
            dash.html.Hr(),
            dcc.Markdown(f"[Wikidata page of {taxon_name}](https://www.wikidata.org/entity/Q{wid})"),
        ]),
        dbc.Row([
            dbc.Col([dbc.Alert(warning, color="primary")]),
            dbc.Col([dbc.Button("Download SMILES", id="btn-download")]),
        ]),
        dcc.Download(id="download"),
        dbc.Row([
            dbc.Pagination(
                id="pagination",
                max_value=math.ceil(nb_matches / PAGE_SIZE),
                fully_expanded=False,
                size="sm",
            ),
        ]),
        dbc.Spinner(id="loading-structures-tsv", children=[dbc.Row(id="cards")]),
    ])
