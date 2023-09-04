from typing import Annotated, List

from fastapi import FastAPI, Query
from fastapi.middleware.wsgi import WSGIMiddleware
from fastapi_versioning import VersionedFastAPI, version
from pydantic import BaseModel
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


class StructureCountResult(BaseModel):
    structure_count: int
    description: str

class StructureResult(BaseModel):
    structure_id: List[int]
    structure_smiles: List[str]
    description: str

class TaxonCountResult(BaseModel):
    taxon_count: int
    description: str

class TaxonResult(BaseModel):
    taxon_id: List[int]
    taxon_name: List[str]
    description: str


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
) -> StructureCountResult:
    results = dm.num_compounds()
    desc = "Number of structures matching the query"
    if q:
        results.update({"q": q})
    return StructureCountResult(structure_count=results, description=desc)


@app.get("/structures/{structure_id}")
@version(1, 0)
async def read_structure(structure_id: str, q: str | None = None, short: bool = False) -> StructureResult:
    ## COMMENT (AR): Make it work for SMILES, InChIKeys, InChIs, names?
    results =  list(dm.compound_search(structure_id))
    ## COMMENT (AR): Throwing out score for now, quite dirty
    results_filtered = [id for id, score in results if score == 1]
    # if q:
    #     results.update({"q": q})
    # if not short:
    #     results.update(
    #         {"description": "This is an amazing item that has a long description"}
    #     )
    structure_smileses = dm.get_compound_smiles_from_list_of_wid(results_filtered)
    desc = "Structures matching the query"
    return StructureResult(structure_id=results_filtered, structure_smiles=structure_smileses, description=desc)

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
) -> TaxonCountResult:
    results = dm.num_taxa()
    desc = "Number of taxa matching the query"
    if q:
        results.update({"q": q})
    return TaxonCountResult(taxon_count=results, description=desc)


@app.get("/taxa/{taxon_id}")
@version(1, 0)
async def read_taxon(taxon_id: str, q: str | None = None, short: bool = False) -> TaxonResult:
    results =  list(dm.get_taxa_with_name_containing(taxon_id))
    # if q:
    #     results.update({"q": q})
    # if not short:
    #     results.update(
    #         {"description": "This is an amazing item that has a long description"}
    #     )
    ## COMMENT (AR): Have a list variant as for the structures?
    ## COMMENT (AR): Make it work with children taxa
    taxon_names = []
    for wid in results:
        taxon_name = dm.get_taxon_name_from_wid(wid)
        if taxon_name is not None:
            taxon_names.append(taxon_name)
    desc = "Taxa matching the query"
    return TaxonResult(taxon_id=results, taxon_name=taxon_names, description=desc)

# @app.post("/items/")
# async def create_item(item: Item):
#     return item


app = VersionedFastAPI(app, enable_latest=True)

app.mount("", WSGIMiddleware(dashboard1.server))
