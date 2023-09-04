from typing import Annotated

from fastapi import FastAPI, Query
from fastapi.middleware.wsgi import WSGIMiddleware
from fastapi_versioning import VersionedFastAPI, version
from dash_app import app as dashboard1
from model import DataModel

description = """
LOTUS API helps you do awesome stuff. 🚀

## Items

### Structures

You can **read structures**.

### Taxa

You can **read taxa**.

## Users

You will be able to:

* **Create users** (_not implemented_).
* **Read users** (_not implemented_).
"""

# from pydantic import BaseModel


# class Item(BaseModel):
#     name: str
#     description: str | None = None
#     price: float
#     tax: float | None = None

dm = DataModel()

app = FastAPI(
    title="LOTUS FastAPI",
    description=description,
    summary="An awesome way to access natural products related data.",
    # contact={
    #     "name": "Deadpoolio the Amazing",
    #     "url": "http://x-force.example.com/contact/",
    #     "email": "dp@x-force.example.com",
    # },
    # license_info={
    #     "name": "Apache 2.0",
    #     "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    # },
        )

@app.get("/structures/")
@version(1, 0)
async def read_structures(
    q: Annotated[
        str | None,
        Query(
            alias="item-query",
            title="Query string",
            description="Query string for the items to search in the database that have a good match",
            min_length=3,
            max_length=50,
            pattern="^fixedquery$"
            ),
    ] = None
):
    results = dm.num_compounds()
    if q:
        results.update({"q": q})
    return results


@app.get("/structures/{structure_id}")
@version(1, 0)
async def read_structure(structure_id: str, q: str | None = None, short: bool = False):
    item = {"structure_id": structure_id}
    if q:
        item.update({"q": q})
    if not short:
        item.update(
            {"description": "This is an amazing item that has a long description"}
        )
    return item

@app.get("/taxa/")
@version(1, 0)
async def read_taxa(
    q: Annotated[
        str | None,
        Query(
            alias="item-query",
            title="Query string",
            description="Query string for the items to search in the database that have a good match",
            min_length=3,
            max_length=50,
            pattern="^fixedquery$"
            ),
    ] = None
):
    results = dm.num_taxa()
    if q:
        results.update({"q": q})
    return results


@app.get("/taxa/{taxon_id}")
@version(1, 0)
async def read_taxon(taxon_id: str, q: str | None = None, short: bool = False):
    item = {"taxon_id": taxon_id}
    if q:
        item.update({"q": q})
    if not short:
        item.update(
            {"description": "This is an amazing item that has a long description"}
        )
    return item

# @app.post("/items/")
# async def create_item(item: Item):
#     return item


app = VersionedFastAPI(app, enable_latest=True)

app.mount("", WSGIMiddleware(dashboard1.server))
