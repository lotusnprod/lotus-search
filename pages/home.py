import dash
import dash_bootstrap_components as dbc
from dash import html

dash.register_page(__name__, path='/', order=1)



def layout():
    return dbc.Container([
        html.Img(src=dash.get_asset_url('lotus.png')),
        dash.dcc.Markdown(open("home.md", "r").read())
    ])
