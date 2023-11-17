import dash
import dash_bootstrap_components as dbc
import plotly
from dash import Dash
from flask import Flask

from model import DataModel

plotly.io.json.config.default_engine = "orjson"

dash._dash_renderer._set_react_version("18.2.0")

data = DataModel()

server = Flask(__name__)
app = Dash(
    __name__,
    server=server,
    external_stylesheets=[dbc.themes.FLATLY],
    use_pages=True,
    suppress_callback_exceptions=True,
)

navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink(page["name"], href=page["path"], active="exact"))
        for page in dash.page_registry.values()
        if page["order"] >= 0
    ],
    brand="Lotus",
    brand_href="/",
    color="primary",
    dark=True,
)

app.layout = dbc.Container([navbar, dash.page_container])

if __name__ == "__main__":
    app.run_server(debug=True, dev_tools_hot_reload=False, use_reloader=False)
