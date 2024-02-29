import dash_bootstrap_components as dbc

import dash
from dash import Input, Output, callback, dash_table, dcc
from model.model import DataModel

dash.register_page(
    __name__, name="Taxon search", top_nav=True, path="/taxa/search", order=1
)

data = DataModel()


@callback(
    Output("taxon-list", "data"),
    Input("input-on-submit", "value"),
)
def update_table(value):
    global data
    if value == "" or len(value) < 3:
        return None

    output = []
    for match in data.get_taxa_with_name_containing(value):
        name = data.get_taxon_object_from_tid(match)
        matching_structures = data.get_structures_of_taxon(match)
        safe_name = name.replace("[", "").replace("]", "")
        output.append(
            {
                "Taxon": f"[{safe_name}](/taxon/{match})",
                "Rank": data.get_ranks_string(match),
                "structures": len(matching_structures),
            }
        )

    return sorted(output, key=lambda x: x["Taxon"])


def layout(name: str = ""):
    return dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Alert(
                        dcc.Markdown(
                            """
            Type a part of the name of the organism you are looking for, for example *ana lutea* will match **Gentiana lutea**.
            
            If you do not find your organism name, use the [Taxon resolver](/taxon_resolver) as it may have a new accepted name"""
                        ),
                        color="success",
                    )
                ]
            ),
            dbc.Row([dbc.Col([dbc.Label("Search for taxon name containing:")])]),
            dbc.Row(
                [dbc.Col([dbc.Input(id="input-on-submit", type="text", value=name)])]
            ),
            dbc.Row([dbc.Col([dbc.Label(" ")])]),
            dbc.Row(
                [
                    dash_table.DataTable(
                        data=None,
                        page_size=15,
                        id="taxon-list",
                        sort_action="native",
                        filter_action="native",
                        columns=[
                            {
                                "name": "Taxon",
                                "id": "Taxon",
                                "type": "text",
                                "presentation": "markdown",
                            },
                            {"name": "Rank ", "id": "Rank"},
                            {"name": "structures", "id": "structures"},
                        ],
                    )
                ]
            ),
        ]
    )
