import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dash_table, dcc

from model import DataModel

dash.register_page(__name__, name="Taxon search", top_nav=True, path="/taxon/search", order=1)


data = DataModel()


def layout():
    return dbc.Container([
        dbc.Input(id='input-on-submit', type='text', value=''),
        dash_table.DataTable(data=None, page_size=15, id='taxon-list',
                             sort_action='native', filter_action='native',
                             columns=[
                                 {'name': 'Taxon', 'id': 'Taxon', 'type': 'text', 'presentation': 'markdown'},
                                 {'name': 'Compounds', 'id': 'Compounds'},
                             ],
                             )
    ])


@callback(
    Output('taxon-list', 'data'),
    Input('input-on-submit', 'value'),
)
def update_table(value):
    global data
    if value == "" or len(value) < 3:
        return None

    output = []
    for match in data.get_taxa_with_name_containing(value):
        name = data.get_taxon_name_from_wid(match)
        matching_compounds = data.get_compounds_of_taxon(match)
        safe_name = name.replace("[", "").replace("]", "")
        output.append(
            {
                "Taxon": f"[{safe_name}](/taxon?wid={match})",
                "Compounds": len(matching_compounds),
            }
        )

    return sorted(output, key=lambda x: x["Taxon"])
