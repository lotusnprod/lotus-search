from urllib.parse import quote

import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

from dash_common import get_svg_of_wid
from model import DataModel

dm = DataModel()


def title(wid=None):
    if wid is not None:
        return f"LOTUS - Q{wid}"

    return "LOTUS"


dash.register_page(
    __name__,
    name="Compound information",
    top_nav=True,
    order=-1,
    path_template="/structure/<wid>",
    title=title,
)


def layout(wid: int):
    if wid is None:
        return dbc.Container([])
    try:
        wid = int(wid)
    except ValueError:
        return dbc.Container([])

    img = get_svg_of_wid(wid)
    img_data = f"data:image/svg+xml,{quote(img)}"

    name_id_list = []
    for t in dm.get_taxa_containing_compound(wid):
        name = dm.get_taxon_name_from_wid(t)
        name_id_list.append([name, t])
    name_id_list = sorted(name_id_list, key=lambda x: x[0])
    table = [{"Taxon": f"[{x[0]}](/taxon/{x[1]})"} for x in name_id_list]
    n_tax = len(name_id_list)
    warning = f"Found in {n_tax} {'taxa' if n_tax > 1 else 'taxon'}"
    return dbc.Container(
        [
            dbc.Row(
                [
                    dash.html.H1(f"Q{wid}"),
                    dash.html.Hr(),
                    dcc.Markdown(
                        f"[Wikidata page of Q{wid}](https://www.wikidata.org/entity/Q{wid})"
                    ),
                    dbc.Row(
                        [
                            dbc.Col([html.Img(src=img_data)]),
                        ]
                    ),
                    dbc.Row([dbc.Alert(warning, color="primary")]),
                    dbc.Row(
                        [
                            dash_table.DataTable(
                                data=table,
                                page_size=15,
                                id="taxon-list-compound",
                                sort_action="native",
                                filter_action="native",
                                columns=[
                                    {
                                        "name": "Taxon",
                                        "id": "Taxon",
                                        "type": "text",
                                        "presentation": "markdown",
                                    },
                                ],
                            )
                        ]
                    ),
                ]
            ),
        ]
    )
