import math

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dcc
from rdkit import Chem

import plotly_dash_ketcher
from config import PAGE_SIZE
from dash_common import generate_compounds_cards
from model import DataModel
from processing_common import fingerprint, standardize

dash.register_page(__name__, name="Structure search", top_nav=True, path="/structures/search", order=2)

dm = DataModel()


def search(fp: bytes) -> list[tuple[int, float]]:
    return dm.compound_search(fp)


def ss_search(fp: bytes, mol, chirality: bool) -> list[tuple[int, float]]:
    return dm.compound_search_substructure(fp, mol, chirality)


def get_matching_ids(query: str,
                     ss_mode: bool, chirality: bool,
                     level: float) -> list[tuple[int, float]]:
    if query:
        try:
            # Sometimes ketcher gives really invalid smiles like with theobromine
            mol = standardize(Chem.MolFromSmiles(query))

            fp = fingerprint(mol)

            if ss_mode:
                scores = ss_search(fp, mol, chirality)
            else:
                scores = search(fp)
                scores = [score for score in scores if score[1] >= level]

            scores_sorted = sorted(scores, reverse=True, key=lambda x: x[1])

            return scores_sorted
        except Exception:
            pass
    return []


@callback(
    Output('chirality-cb', 'disabled'),
    Input('substructure-cb', 'value')
)
def disable_chirality(value):
    return not value


@callback(
    Output('similarity-slider', 'disabled'),
    Input('substructure-cb', 'value')
)
def disable_similarity(value):
    return value


@callback(
    Output("download-compounds", "data"),
    [Input("btn-download-compounds", "n_clicks"),
     Input("structure-search-data", "data")],
    prevent_initial_call=True,
)
def download_compounds(n_clicks, data):
    if dash.ctx.triggered_id == "btn-download-compounds":
        filename = f"compounds.tsv"
        return dict(content=dm.compound_get_tsv_from_scores(data["matching_compounds"], data["scores"]),
                    filename=filename)


@callback(
    Output("structure-search-data", "data"),
    Output("structure-search-results", "children"),
    Output("pagination-compound-search", "max_value"),
    Output("warning-compound-search", "children"),
    [Input("ketcher", "molecule"),
     Input("substructure-cb", "value"),
     Input("chirality-cb", "value"),
     Input("similarity-slider", "value"),
     Input("pagination-compound-search", "active_page")]
)
def search_compound_cards(molecule: str,
                          ss_mode: bool,
                          chirality: bool,
                          similarity: float,
                          active_page: int):
    data = {}
    scores = get_matching_ids(molecule, ss_mode, chirality, similarity)
    n_scores = len(scores)
    if n_scores > 0:
        warning = f"Found {n_scores} matching compounds"
    else:
        warning = f"No matching compounds found"
    data["matching_compounds"] = [score[0] for score in scores]
    data["scores"] = [score[1] for score in scores]
    return (data,
            generate_compounds_cards(active_page, data),
            math.ceil(len(scores) / PAGE_SIZE),
            warning)


def layout():
    return dbc.Container([
        dcc.Store(id='structure-search-data', data={
            "matching_compounds": [],
            "scores": []
        }),
        plotly_dash_ketcher.PlotlyDashKetcher(id="ketcher", buttonLabel="Search"),
        dbc.Row([
            dbc.Col([
                dbc.Checkbox(id='substructure-cb', value=False, label='Substructure search'),
            ]),
            dbc.Col([
                dbc.Checkbox(id='chirality-cb', value=False, disabled=True,
                             label='Chirality search (has defects)'),
            ])
        ]),
        dbc.Row([
            dbc.Label("Similarity level"),
            dbc.Col([
                dcc.Slider(id="similarity-slider",
                           min=0, max=1, step=0.01, marks=None, value=0.8)
            ])
        ]),
        dbc.Row([
            dbc.Col([dbc.Alert(id="warning-compound-search", color="primary")]),
            dbc.Col(
                dbc.Button("Download SMILES", id="btn-download-compounds")
            )
        ]),
        dbc.Spinner(id="loading-compounds-tsv",
                    children=[dcc.Download(id="download-compounds")],
                    ),
        dbc.Row([
            dbc.Col(dbc.Pagination(id="pagination-compound-search", max_value=1, fully_expanded=False,
                                   size="lg", first_last=True, previous_next=True)),
        ]),
        dbc.Spinner(id="loading-compounds-search",
                    children=[dbc.Row(id='structure-search-results')])
    ])