import json
import time
import uuid

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback, dcc, html

import plotly_dash_ketcher
from model import DataModel

dm = DataModel()

dash.register_page(__name__, name="Something is missing", path='/contribute', order=90)


@callback(
    Output("organism_contribution", "disabled"),
    Output("organism_part_row", "style"),
    Input("ketcher-contrib", "molecule")
)
def enable_organism(molecule):
    disabled = molecule is None or molecule == ""
    return disabled, {"display": "block" if not disabled else "none"}


@callback(
    Output("doi", "disabled"),
    Output("search_organism_contribution", "disabled"),
    Output("doi_part_row", "style"),
    Input("organism_contribution", "value")
)
def enable_doi(value):
    disabled = value is None or value == "" or len(value) < 3
    return disabled, disabled, {"display": "block" if not disabled else "none"}


@callback(
    Output("submit_contribution_single", "disabled"),
    Output("submit_part_row", "style"),
    Input("doi", "value")
)
def enable_submit(value):
    disabled = value is None or value == ""
    return disabled, {"display": "block" if not disabled else "none"}


@callback(
    Output("matching_organisms_contribution_select", "options"),
    Input("organism_contribution", "value"),
    Input("search_organism_contribution", "n_clicks")
)
def search_organism(name, n_clicks):
    if dash.ctx.triggered_id == "search_organism_contribution":
        results = dm.get_taxa_with_name_containing(name)
        names = sorted([dm.get_taxon_name_from_wid(match) for match in results])
        return [{"label": name, "value": name} for name in names]
    return []


tab1_content = dbc.Card(
    dbc.CardBody(
        [
            dbc.Label("Draw the structure and press 'Add'"),
            plotly_dash_ketcher.PlotlyDashKetcher(id="ketcher-contrib", buttonLabel="Add"),
            dbc.Form(id="organism_part_row", children=[
                dbc.Label("Organism name (mandatory)", html_for="organism_contribution"),
                dbc.Input(id="organism_contribution", placeholder="Enter the organism name", disabled=True),
                dbc.Button(id="search_organism_contribution", children=["Search organism"], color="primary",
                           disabled=True),
                dbc.Select(id="matching_organisms_contribution_select",
                           options=[],
                           placeholder="Click on Search organism and select organism from this list (optional)")]),
            html.Div(id="doi_part_row", children=[
                dbc.Label("Publication reference (mandatory)", html_for="doi"),
                dbc.Input(id="doi", placeholder="Enter the publication reference (DOI is prefered)", disabled=True)
            ]),
            html.Div(id="submit_part_row", children=[dbc.Button(id="submit_contribution_single", children=["Submit"],
                                                                color="success", disabled=True),
                                                     dbc.Row(id="result_contribution_single")])
        ]
    ),
    className="mt-3",
)


@callback(
    Output("result_contribution_single", "children"),
    Input("submit_contribution_single", "n_clicks"),
    Input("email", "value"),
    Input("ketcher-contrib", "molecule"),
    Input("organism_contribution", "value"),
    Input("matching_organisms_contribution_select", "value"),
    Input("doi", "value")
)
def submit(n_clicks, email, structure, organism, selected_organism, reference):
    if dash.ctx.triggered_id == "submit_contribution_single":
        if email is None:
            email = ""
        if structure is None:
            structure = ""
        if organism is None:
            organism = ""
        if reference is None:
            reference = ""
        if selected_organism is None:
            selected_organism = ""

        if len(structure) > 8192:
            return dbc.Alert("Your molecule is too big, we do not handle molecules this size yet.", color="danger")
        elif len(organism) > 100 or len(selected_organism) > 100:
            return dbc.Alert("Your organism name is too long, please use a shorter name.", color="danger")
        elif len(reference) > 512:
            return dbc.Alert("Your reference is way too long, please use a doi or a short reference.", color="danger")
        elif len(email) > 64 or ("@" not in email and email != ""):
            return dbc.Alert(f"Your email {email} is invalid, please use a valid email address or none at all",
                             color="danger")

        js = json.dumps({
            "molecule": structure,
            "organism": organism,
            "selected_organism": selected_organism,
            "reference": reference,
            "email": email
        })
        with open(f"contributions/{time.time()}_{uuid.uuid4()}.json", "w") as f:
            f.write(js)
        return dbc.Alert("Thank you for your contribution, we will get back to you as soon as possible.",
                         color="success")
    return None


@callback(Output('output-data-upload', 'children'),
          Input("email", "value"),
          Input('upload-data', 'contents'),
          State('upload-data', 'filename'))
def update_output(email, content, name):
    if content is not None:
        if len(content) > 1024 * 1024:
            return dbc.Alert(f"Sorry your file {name} is too big, we do not handle files this size yet.",
                             color="danger")
        if not name.endswith(".csv") or not name.endswith(".tsv"):
            return dbc.Alert(f"Sorry your file {name} is not a CSV or TSV file.", color="danger")
        basename = f"contributions/{time.time()}_{uuid.uuid4()}"
        with open(f"{basename}.xsv", "w") as f:
            f.write(content)
        if email != "" and email is not None and "@" in email and len(email) < 64:
            with open(f"{basename}.txt", "w") as f:
                f.write(email)
        return dbc.Alert(
            f"Thank you for your contribution, we received {name} successfully, we will get back to you as soon as possible.",
            color="success")
    return None


tab2_content = dbc.Card(
    dbc.CardBody(
        [
            dcc.Markdown("""You are even more awesome!

So if you have way too many compounds to add, you can send us directly a TSV or CSV file with the following three columns:

**smiles** , **organism** , **reference**

The reference should ideally be a DOI. And the organism should ideally be the latest accepted name.

We will only be able to add from articles we have access to (we have accesses from several universities if necessary).

Priority will be given to articles that are open access.

Please make sure that it is formatted properly and that the values are quoted (even if they don't contain commas)""",
                         className="card-text"),
            dcc.Upload(
                id='upload-data',
                children=dbc.Button([
                    'Drag and Drop or Select Files'
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'solid',
                    'borderColor': '#ccc',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
            dbc.Row(id="output-data-upload")
        ]
    ),
    className="mt-3",
)


def layout():
    return dbc.Container([
        html.H1("Awesome, lets add it"),
        dbc.Row([html.Div(
            "That means you will be able to help us improve the data quality of LOTUS. We trust that you have checked that it is not already present in the database and that your molecule is drawn correctly.")]),
        html.Hr(),
        dbc.Form([html.Div([
            dbc.Label("Your email (optional)", html_for="email"),
            dbc.Input(id="email", type="email", placeholder="Enter email (optional)"),
            dbc.FormText(
                "We will contact you if we have questions, or to let you know it is added. Your contribution may also be prioritized.",
                color="secondary",
            )])
        ]),
        dbc.Row([dbc.Tabs(
            [
                dbc.Tab(tab1_content, label="I have a single compound"),
                dbc.Tab(tab2_content, label="I have multiple compounds"),
            ]
        )])
    ])
