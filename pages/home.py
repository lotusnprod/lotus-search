import dash_bootstrap_components as dbc

import dash
from dash import html

dash.register_page(__name__, path="/", order=1)


def layout():
    return dbc.Container([
        html.Img(
            src="https://upload.wikimedia.org/wikipedia/commons/6/64/Lotus_initiative_logo.svg",
        ),
        dash.dcc.Markdown(open("home.md").read()),
    ])
