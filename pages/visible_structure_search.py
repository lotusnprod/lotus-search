import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dash_table
import plotly_dash_ketcher
from model import DataModel

dash.register_page(__name__, name="Structure search", top_nav=True, path="/structures/search", order=2)


data = DataModel()

@callback(
    Output('test-smiles', 'value'),
    Input('component', 'molecule'),
)
def update_smiles(molecule):
    return molecule

@callback(
    Output('chirality-cb', 'disabled'),
    Input('substructure-cb', 'value')
)
def disable_chirality(value):
    return not value

def layout():
    return dbc.Container([
        plotly_dash_ketcher.PlotlyDashKetcher(id='component'),
        dbc.Input(id='test-smiles', type='text', value=''),
        dbc.Row([
            dbc.Col([
                dbc.Checkbox(id='substructure-cb', value=False, label='Substructure search'),
            ]),
            dbc.Col([
                dbc.Checkbox(id='chirality-cb', value=False, disabled=True,
                             label='Chirality search (has defects)'),
            ])
        ])
    ])

