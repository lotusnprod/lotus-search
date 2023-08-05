import dash
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/', order=1)



def layout():
    return dbc.Container([
        dash.dcc.Markdown(open("home.md", "r").read())
    ])
