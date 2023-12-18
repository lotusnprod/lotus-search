import math

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dcc

import plotly_dash_ketcher
from dash_config import PAGE_SIZE
from dash_common import generate_structures_cards
from model import DataModel

dash.register_page(
    __name__, name="Structure search", top_nav=True, path="/structures/search", order=2
)

dm = DataModel()


def get_matching_ids(
    query: str, ss_mode: bool, chirality: bool, level: float
) -> list[tuple[int, float]]:
    if query:
        try:
            # Sometimes ketcher gives really invalid smiles like with theobromine
            # TODO move all that in the model

            if ss_mode:
                scores = dm.structure_search_substructure(query, chirality)
            else:
                scores = dm.structure_search(query)
                scores = [score for score in scores if score[1] >= level]

            scores_sorted = sorted(scores, reverse=True, key=lambda x: x[1])

            return scores_sorted
        except Exception:
            pass
    return []


@callback(Output("chirality-cb", "style"), Input("substructure-cb", "value"))
def disable_chirality(value):
    return {"display": "block" if value else "none"}


@callback(Output("similarity-row", "style"), Input("substructure-cb", "value"))
def disable_similarity(value):
    return {"display": "none" if value else "block"}


@callback(
    Output("download-structures", "data"),
    [
        Input("btn-download-structures", "n_clicks"),
        Input("structure-search-data", "data"),
    ],
    prevent_initial_call=True,
)
def download_structures(n_clicks, data):
    if dash.ctx.triggered_id == "btn-download-structures":
        filename = "structures.tsv"
        return dict(
            content=dm.structure_get_tsv_from_scores(
                data["matching_structures"], data["scores"]
            ),
            filename=filename,
        )


@callback(
    Output("structure-search-data", "data"),
    Output("structure-search-results", "children"),
    Output("pagination-structure-search", "max_value"),
    Output("pagination-structure-search", "style"),
    Output("btn-download-structures", "style"),
    Output("warning-structure-search", "children"),
    [
        Input("ketcher", "molecule"),
        Input("substructure-cb", "value"),
        Input("chirality-cb", "value"),
        Input("similarity-slider", "value"),
        Input("pagination-structure-search", "active_page"),
    ],
)
def search_structure_cards(
    molecule: str, ss_mode: bool, chirality: bool, similarity: float, active_page: int
):
    data = {}
    scores = get_matching_ids(molecule, ss_mode, chirality, similarity)
    n_scores = len(scores)
    if n_scores > 0:
        warning = f"Found {n_scores} matching structures"
    else:
        warning = "No matching structures found"
    data["matching_structures"] = [score[0] for score in scores]
    data["scores"] = [score[1] for score in scores]
    style_pagination = {"display": "none" if n_scores == 0 else "flex"}
    style_btn = {"display": "none" if n_scores == 0 else "block"}
    return (
        data,
        generate_structures_cards(active_page, data, molecule),
        math.ceil(len(scores) / PAGE_SIZE),
        style_pagination,
        style_btn,
        warning,
    )


def layout():
    return dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Alert(
                        dcc.Markdown(
                            """
Draw a molecule (or substructure) below.
You can also simply copy/paste a SMILES (such as `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` for caffeine).
                    """
                        ),
                        color="success",
                    )
                ]
            ),
            dcc.Store(
                id="structure-search-data",
                data={"matching_structures": [], "scores": []},
            ),
            plotly_dash_ketcher.PlotlyDashKetcher(id="ketcher", buttonLabel="Search"),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Checkbox(
                                id="substructure-cb",
                                value=False,
                                label="Substructure search",
                            ),
                        ]
                    ),
                    dbc.Col(
                        [
                            dbc.Checkbox(
                                id="chirality-cb",
                                value=False,
                                label="Chirality search (has defects)",
                            ),
                        ]
                    ),
                ]
            ),
            dbc.Row(
                id="similarity-row",
                children=[
                    dbc.Label("Similarity level"),
                    dbc.Col(
                        [dcc.Slider(id="similarity-slider", min=0, max=1, value=0.8)]
                    ),
                ],
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [dbc.Alert(id="warning-structure-search", color="primary")]
                    ),
                    dbc.Col(
                        dbc.Button("Download SMILES", id="btn-download-structures")
                    ),
                ]
            ),
            dbc.Spinner(
                id="loading-structures-tsv",
                children=[dcc.Download(id="download-structures")],
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Pagination(
                            id="pagination-structure-search",
                            max_value=1,
                            fully_expanded=False,
                            size="lg",
                            first_last=True,
                            previous_next=True,
                        )
                    ),
                ]
            ),
            dbc.Spinner(
                id="loading-structures-search",
                children=[dbc.Row(id="structure-search-results")],
            ),
        ]
    )
