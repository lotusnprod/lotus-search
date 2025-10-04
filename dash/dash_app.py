import dash_bootstrap_components as dbc
import plotly
from flask import Flask

import dash
from dash import Dash, html
from dash.data_provider import get_data_model

plotly.io.json.config.default_engine = "orjson"

dash._dash_renderer._set_react_version("18.2.0")

data = get_data_model()

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
    class_name="mb-3",
)

footer = dbc.Container(
    [
        html.Hr(),
        html.Small(
            "LOTUS Search – Experimental. Data © LOTUS / Wikidata. ",
            className="text-muted",
        ),
    ],
    fluid=True,
    className="mt-4 mb-2",
)

app.layout = dbc.Container([navbar, dash.page_container, footer], fluid=True)

if __name__ == "__main__":
    app.run_server(debug=True, dev_tools_hot_reload=False, use_reloader=False)
