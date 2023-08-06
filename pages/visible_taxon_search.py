import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dash_table

from model import DataModel

dash.register_page(__name__, name="Taxon search", top_nav=True, path="/taxa/search", order=1)

data = DataModel()


@callback(
    Output('taxon-list', 'data'),
    Input('input-on-submit', 'value'),
)
def update_table(value):
    global data
    if value == "" or len(value) < 3:
        return None

    output = []
    for match in data.get_taxa_with_name_containing(value):
        name = data.get_taxon_name_from_wid(match)
        matching_compounds = data.get_compounds_of_taxon(match)
        safe_name = name.replace("[", "").replace("]", "")
        output.append(
            {
                "Taxon": f"[{safe_name}](/taxon/{match})",
                "Compounds": len(matching_compounds),
            }
        )

    return sorted(output, key=lambda x: x["Taxon"])


def layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([dbc.Label("Search for taxon name containing:")])]),
        dbc.Row([
            dbc.Col([dbc.Input(id='input-on-submit', type='text', value='')
                     ])])
        ,
        dbc.Row([dbc.Col([dbc.Label(" ")])]),
        dbc.Row([
        dash_table.DataTable(data=None, page_size=15, id='taxon-list',
                             sort_action='native', filter_action='native',
                             columns=[
                                 {'name': 'Taxon', 'id': 'Taxon', 'type': 'text', 'presentation': 'markdown'},
                                 {'name': 'Compounds', 'id': 'Compounds'},
                             ],
                             )
        ])
    ])
