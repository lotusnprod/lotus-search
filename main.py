from typing import Annotated, Dict, List

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

class CoupleDict(BaseModel):
    structure_id: List[int]
    structure_smiles: List[str]
    taxon: Dict

class CoupleResult(BaseModel):
    results: CoupleDict
    count: int
    description: str

class StructureDict(BaseModel):
    structure_id: List[int]
    structure_smiles: List[str]

class StructureResult(BaseModel):
    results: StructureDict
    count: int
    description: str
    
class TaxonDict(BaseModel):
    taxon: Dict

class TaxonResult(BaseModel):
    results: TaxonDict
    count: int
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

@app.post("/couples/")
@version(1, 0)
async def create_couples(
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
) -> CoupleResult:
    desc = "Couples matching the query"
    ## TODO
    return CoupleResult(results=CoupleDict(structure_id=results, structure_smiles=smiles), description=desc, count = len(results))


@app.post("/couples/{couples_query}") # no example for now
@version(1, 0)
async def create_couples(couple_query: str, q: str | None = None) -> CoupleResult:
    ## TODO
    return CoupleResult(results=CoupleDict(), description=desc, count = len(ids_filtered))

@app.post("/structures/")
@version(1, 0)
async def create_structures(
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
) -> StructureResult:
    desc = "Structures matching the query"
    results = dm.get_compounds()
    if q:
        results.update({"q": q})
    ## For dev
    results = results[:500]
    smiles = dm.get_compound_smiles_from_list_of_wid(results)

    return StructureResult(results=StructureDict(structure_id=results, structure_smiles=smiles), description=desc, count = len(results))


@app.post("/structures/{structure_query}") # call the id WID? -> No, no ID is intended here, renamed it (example: http://127.0.0.1:8000/v1_0/structures/CC1C=C(C(=O)C2(C1CC3C4(C2C(=O)C(=C(C4CC(=O)O3)C)OC)C)C)OC)
@version(1, 0)
async def create_structure(structure_query: str, q: str | None = None) -> StructureResult:
    desc = "Structures matching the query"
    ## COMMENT (AR): Make it work for SMILES, InChIKeys, InChIs, names?
    ids =  list(dm.compound_search(structure_query))
    ## COMMENT (AR): Throwing out score for now, quite dirty
    ids_filtered = [id for id, score in ids if score == 1]
    # if q:
    #     results.update({"q": q})
    structure_smileses = dm.get_compound_smiles_from_list_of_wid(ids_filtered)

    return StructureResult(results=StructureDict(structure_id=ids_filtered, structure_smiles=structure_smileses), description=desc, count = len(ids_filtered))

@app.post("/taxa/")
@version(1, 0)
async def create_taxa(
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
) -> TaxonResult:
    desc = "Taxa matching the query"
    results = dm.get_taxa()
    if q:
        results.update({"q": q})
    ## For dev
    results = dict(list(results.items())[:500])

    return TaxonResult(results=TaxonDict(taxon=results), description=desc, count = len(results))


@app.post("/taxa/{taxon_query}")  # call the id WID? -> No, no ID is intended here, renamed it (example: http://127.0.0.1:8000/v1_0/taxa/gentiano)
@version(1, 0)
async def create_taxon(taxon_query: str, q: str | None = None) -> TaxonResult:
    desc = "Taxa matching the query"
    ids = list(dm.get_taxa_with_name_containing(taxon_query))
    # if q:
    #     results.update({"q": q})
    taxon_names = dm.get_taxon_name_from_list_of_wid(ids)
    results = dict(zip(ids, taxon_names))

    return TaxonResult(results=TaxonDict(taxon=results), description=desc, count = len(results))


# The API should not be able to create items, only read them.
# We can have another API or authentication later for that
# @app.post("/items/")
# async def create_item(item: Item):
#     return item


app = VersionedFastAPI(app, enable_latest=True)

app.mount("", WSGIMiddleware(dashboard1.server))
