import dash
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/', order=1)



def layout():
    return dbc.Container([
        dbc.Alert("Hello, world!", color="primary"),
    ])
