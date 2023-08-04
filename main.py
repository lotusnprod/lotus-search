import dash
import dash_bootstrap_components as dbc
from dash import Dash

from model import DataModel

data = DataModel()

app = Dash(__name__,
           external_stylesheets=[dbc.themes.CERULEAN],
           use_pages=True,
           suppress_callback_exceptions=True)

navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink(page["name"], href=page["path"], active="exact"))
        for page in dash.page_registry.values()
        if page["order"] >= 0
    ],
    brand="Lotus", brand_href="/", color="primary", dark=True,
)

app.layout = dbc.Container([
    navbar,
    dash.page_container
])
if __name__ == '__main__':
    app.run(debug=True, threaded=True)
