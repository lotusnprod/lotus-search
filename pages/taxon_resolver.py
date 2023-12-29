from typing import Any
from urllib.parse import quote

import dash_bootstrap_components as dbc

import dash
from dash import Input, Output, State, callback, dcc, html
from model.model import DataModel

dash.register_page(
    __name__, name="Taxon resolver", top_nav=True, path="/taxon_resolver", order=70
)

dm = DataModel()


@callback(
    Output("taxon-resolver-list", "children"),
    Output("taxon-resolver-alert", "children"),
    [
        State("input-on-submit-taxon-resolver", "value"),
        Input("submit-button-taxon-resolver", "n_clicks"),
        Input("taxon-resolver-form", "n_submit"),
    ],
)
def search_taxon(query: str, n_clicks: int, n_submit: int) -> Any:
    print("Triggering")
    if query is None or len(query) < 3:
        return [], "Please enter at least 3 characters."

    try:
        data = dm.resolve_taxon(query=query)
    except Exception as e:
        print(f"Big fail while requesting GN {e}")
        return [], dcc.Markdown(f"Error while requesting GN, please retry. ({e})")
    table_header = [
        html.Thead(
            html.Tr(
                [html.Th("Accepted name"), html.Th("Database"), html.Th("Match type")]
            )
        )
    ]
    outputs = []
    try:
        if "names" not in data:
            print("Error while parsing GN data")
            print(data)
            return [], dcc.Markdown("Error while parsing, please retry.")
        if len(data["names"]) == 0:
            return [], dcc.Markdown(f"Sorry no match found for {query}.")
        for n in data["names"]:
            if "results" not in n or len(n["results"]) == 0:
                return [], dcc.Markdown(f"Sorry no match found for {query}.")
            for d in n["results"]:
                match_type = d["matchType"]
                if "currentName" not in d:
                    continue
                name = d["currentCanonicalSimple"]
                full_name = d["currentName"]
                db_name = d["dataSourceTitleShort"]
                if "outlink" in d:
                    db_match = f"[{db_name} {full_name}]({d['outlink']})"
                else:
                    db_match = f"{db_name} {full_name}"

                outputs.append(
                    html.Tr(
                        [
                            html.Td(
                                dcc.Markdown(
                                    f"[{name}](/taxa/search?name={quote(name)})"
                                )
                            ),
                            html.Td(dcc.Markdown(db_match)),
                            html.Td(dcc.Markdown(match_type)),
                        ]
                    )
                )

        return table_header + [html.Tbody(outputs)], ""
    except Exception as e:
        print(f"Error while parsing GN data: {e}")
        print(data)
        return [], dcc.Markdown(f"Error while parsing, please retry. ({e})")


def layout():
    return dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Alert(
                        dcc.Markdown(
                            """
This service is using the [Global Names Resolver API](https://verifier.globalnames.org). Give them a visit [https://globalnames.org]

It supports fuzzy matching, for example *Jentiana lutea* will match *Gentiana lutea*.

It does not support vernacular names (E.g. tomato) so you will need to use the systematic name (*Solanum lycopersicum*)
                    """
                        ),
                        color="success",
                    )
                ]
            ),
            dbc.Form(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Label(
                                        "Search for taxon:",
                                        html_for="input-on-submit-taxon-resolver",
                                    )
                                ]
                            )
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Input(
                                        id="input-on-submit-taxon-resolver",
                                        type="text",
                                        value="",
                                    )
                                ]
                            )
                        ]
                    ),
                ],
                id="taxon-resolver-form",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Button(
                                id="submit-button-taxon-resolver",
                                children="Search",
                                color="primary",
                            )
                        ]
                    )
                ]
            ),
            html.Div(id="taxon-resolver-alert"),
            dbc.Spinner(
                id="loading-taxa-resolved",
                children=dbc.Table(id="taxon-resolver-list", children=[]),
            ),
        ]
    )
