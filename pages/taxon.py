import math
from urllib.parse import quote

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dcc, html

from chemistry_helpers import molecule_svg
from model import DataModel

dm = DataModel()


def title(wid=None):
    if wid is not None:
        try:
            return f"LOTUS - {dm.get_taxon_name_from_wid(int(wid))}"
        except:
            pass
    return "LOTUS"


def tsv(compounds: list[int]) -> str:
    smileses = dm.get_compound_smiles_from_list_of_wid(compounds)
    return "smiles\n" + "\n".join(smileses)


dash.register_page(__name__, name="Taxon information",
                   top_nav=True, order=-1, path_template="/taxon/<wid>",
                   title=title)


@callback(
    Output("cards", "children"),
    [Input("pagination", "active_page"),
     Input("matching-ids", "data")],
)
def compound_cards(active_page: int, data: list[int]):
    cards = []
    if active_page is None:
        active_page = 1

    displayed = data["matching_compounds"][100 * (active_page - 1):100 * active_page]
    data_dl = displayed
    # if len(matching_compounds) >= 100:
    #     if c2.checkbox("I want to download them all and I understand it can be really slow"):
    #         data_dl = matching_compounds
    #
    # c1.download_button(f"Download {len(data_dl)} smiles as TSV", tsv(data_dl),
    #                    f"Compounds from {name.replace('.', '')}.tsv",
    #                    "text/tab-separated-values")

    for idx, j in enumerate(displayed):
        img = molecule_svg(dm.get_compound_smiles_from_wid(j))
        img_data = f"data:image/svg+xml,{quote(img)}"
        taxa_count = dm.get_number_of_taxa_containing_compound(j)
        card = dbc.Card([
            dbc.CardImg(src=img_data, top=True),
            dbc.CardBody(
                [
                    html.H4(f"Q{j}", className="card-title"),
                    html.P(
                        f"Found in {taxa_count} {'taxa' if taxa_count > 1 else 'taxon'}",
                        className="card-text",
                    ),
                    dcc.Markdown(f"[Wikidata page of compound](https://www.wikidata.org/entity/Q{j})"),
                    dbc.Button("Compound page", color="primary", href=f"/molecule/{j}"),
                ]
            ),
        ],
            style={"width": "18rem"}, )
        cards.append(card)

    return [*cards]

@callback(
    Output("download", "data"),
    [Input("btn-download", "n_clicks"),
    Input("matching-ids", "data")],
    prevent_initial_call=True,
)
def func(n_clicks, data):
    filename = f"compounds_of_{data['taxon_name'].replace('.', '')}.tsv"
    return dict(content=tsv(data["matching_compounds"]), filename=filename)


def layout(wid: int):
    if wid is None:
        return dbc.Container([])
    try:
        wid = int(wid)
    except ValueError:
        return dbc.Container([])

    taxon_name = dm.get_taxon_name_from_wid(wid)
    parent_ranks = dm.get_ranks_string(wid)

    taxonomic_info = ""
    if wid in dm.db["taxonomy_direct_parents"]:
        markdown = ""
        parent_taxa = dm.db["taxonomy_direct_parents"][wid]
        tree = []
        for parent in parent_taxa:
            tree.append([parent, 1])
            if parent in dm.db["taxonomy_parents_with_distance"]:
                for relative in dm.db["taxonomy_parents_with_distance"][parent]:
                    distance = dm.db["taxonomy_parents_with_distance"][parent][relative]
                    tree.append([relative, distance])
        tree = sorted(tree, key=lambda x: x[1])
        for parent in tree:
            ranks = dm.get_ranks_string(parent[0])

            tax_name = dm.get_taxon_name_from_wid(parent[0])
            if tax_name is not None:
                markdown += f"[{tax_name}{ranks}](/taxon/{parent[0]}) > "
        taxonomic_info = markdown.strip("> ")

    matching_compounds = dm.get_compounds_of_taxon(wid)
    matching_compounds = list(matching_compounds)
    matching_compounds.sort()
    nb_matches = len(matching_compounds)

    warning = f"Found {nb_matches} compounds"

    return dbc.Container([
        dcc.Store(id='matching-ids', data={"matching_compounds": matching_compounds,
                                           "taxon_name": taxon_name}),
        dbc.Row([
            dash.html.H1(f"{taxon_name}{parent_ranks}"),
            dcc.Markdown(taxonomic_info),
            dash.html.Hr(),
            dcc.Markdown(f"[Wikidata page of {taxon_name}](https://www.wikidata.org/entity/Q{wid})"),
        ]),
        dbc.Row([
            dbc.Col([dbc.Alert(warning, color="primary")]),
            dbc.Col([dbc.Button("Download SMILES", id="btn-download")])
            ]),
        dcc.Download(id="download"),
        dbc.Row([
            dbc.Pagination(id="pagination", max_value=math.ceil(nb_matches / 100), fully_expanded=False, size="sm"),
        ]),
        dbc.Row(id="cards"),
    ])
