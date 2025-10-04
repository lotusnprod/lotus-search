import dash_bootstrap_components as dbc

import dash
from dash import html, dcc

dash.register_page(__name__, path="/", order=1)


def layout():
    return dbc.Container([
        html.Img(
            src="https://upload.wikimedia.org/wikipedia/commons/6/64/Lotus_initiative_logo.svg",
        ),
        dash.dcc.Markdown(open("home.md").read()),
        html.Hr(),
        html.H4("Quick links"),
        dbc.ListGroup([
            dbc.ListGroupItem(dcc.Link("Structure search", href="/structures/search")),
            dbc.ListGroupItem(dcc.Link("Taxon search", href="/taxa/search")),
            # dbc.ListGroupItem(dcc.Link("API usage", href="/api")),
            dbc.ListGroupItem(dcc.Link("About", href="/about")),
            dbc.ListGroupItem(dcc.Link("Contribute", href="/contribute")),
        ]),
    ])
