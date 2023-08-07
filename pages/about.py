import dash
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/about', order=99)


def layout():
    return dbc.Container([
        dash.dcc.Markdown(open("README.md", "r").read())
    ])
