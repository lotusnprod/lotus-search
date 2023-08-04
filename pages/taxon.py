import dash
import dash_bootstrap_components as dbc

from model import DataModel

dash.register_page(__name__, name="Taxon information", top_nav=True, order=-1, path="/taxon")


data = DataModel()


def layout(wid: int):
    return dbc.Container([
        dbc.Alert(f"Hello, world looking at {wid}!", color="primary"),
    ])

