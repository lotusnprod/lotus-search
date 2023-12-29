from typing import Any
from urllib.parse import quote

import dash_bootstrap_components as dbc
from dash_config import PAGE_SIZE
from flask_caching import Cache

from chemistry_helpers import molecule_svg
from dash import dcc, get_app, html
from model.model import DataModel

dm = DataModel()

cache = Cache(
    get_app().server,
    config={
        # try 'filesystem' if you don't want to setup redis
        "CACHE_TYPE": "simple",  # That should be fine everybody has the same
    },
)


@cache.memoize(timeout=3600)
def get_svg_of_wid(j: int, molecule: str | None = None) -> str:
    return molecule_svg(dm.get_structure_smiles_from_sid(j), molecule)


@cache.memoize(timeout=3600)
def get_number_of_taxa_for_structure(j: int) -> int:
    return dm.get_number_of_taxa_containing_structure(j)


def generate_structures_cards(
    active_page: int, data: dict[str, Any], molecule: str | None = None
) -> list[dbc.Card]:
    cards = []
    if active_page is None:
        active_page = 1

    scores_mode = "scores" in data

    displayed = data["matching_structures"][
        PAGE_SIZE * (active_page - 1) : PAGE_SIZE * active_page
    ]

    if scores_mode:
        scores = data["scores"][PAGE_SIZE * (active_page - 1) : PAGE_SIZE * active_page]

    for idx, j in enumerate(displayed):
        extras = []
        if scores_mode:
            extras.append(dbc.Row([dbc.Label(f"Similarity {scores[idx]:.2f}")]))
            extras.append(
                dbc.Row(
                    dbc.Col(
                        dbc.Progress(
                            key=f"score_{idx}", value=int(scores[idx] * 100), max=100.0
                        )
                    )
                )
            )

        img = get_svg_of_wid(j, molecule)
        img_data = f"data:image/svg+xml,{quote(img)}"
        taxa_count = get_number_of_taxa_for_structure(j)
        card = dbc.Card(
            [
                dbc.CardImg(src=img_data, top=True),
                dbc.CardBody(
                    [
                        html.H4(f"Q{j}", className="card-title"),
                        html.P(
                            f"Found in {taxa_count} {'taxa' if taxa_count > 1 else 'taxon'}",
                            className="card-text",
                        ),
                        dcc.Markdown(
                            f"[Wikidata page of structure](https://www.wikidata.org/entity/Q{j})"
                        ),
                        dbc.Button(
                            "structure page", color="primary", href=f"/structure/{j}"
                        ),
                        *extras,
                    ]
                ),
            ],
            style={"width": "18rem"},
        )
        cards.append(card)

    return cards
