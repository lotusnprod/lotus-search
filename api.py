from typing import Annotated, Dict, List

from fastapi import FastAPI, Query
from fastapi_versioning import VersionedFastAPI, version
from pydantic import BaseModel
from model import DataModel

description = """
LOTUSFast API helps you do awesome stuff. 🚀
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
            alias="query",
            title="Query string",
            description="Query string for the items to search in the database that have a good match"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None
) -> CoupleResult:
    desc = "Couples matching the query"
    ## TODO
    return CoupleResult(results=CoupleDict(structure_id=results, structure_smiles=smiles), description=desc, count = len(results))

@app.post("/structures/")
@version(1, 0)
async def create_structures(
    q: Annotated[
        str | None,
        Query(
            alias="query",
            title="Query string",
            description="Query string for the items to search in the database that have a good match"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None
) -> StructureResult:
    desc = "Structures matching the query"
    results = dm.get_compounds()
    if q:
        ids =  list(dm.compound_search(q))
        ## COMMENT (AR): Throwing out score for now, quite dirty
        results = [id for id, score in ids if score == 1]
    ## For dev
    results = results[:500]
    smiles = dm.get_compound_smiles_from_list_of_wid(results)

    return StructureResult(results=StructureDict(structure_id=results, structure_smiles=smiles), description=desc, count = len(results))

@app.post("/taxa/")
@version(1, 0)
async def create_taxa(
    q: Annotated[
        str | None,
        Query(
            alias="query",
            title="Query string",
            description="Query string for the items to search in the database that have a good match"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None
) -> TaxonResult:
    desc = "Taxa matching the query"
    results = dm.get_taxa()
    if q:
        ids = list(dm.get_taxa_with_name_containing(q))
        taxon_names = dm.get_taxon_name_from_list_of_wid(ids)
        results = dict(zip(ids, taxon_names))

    ## For dev
    results = dict(list(results.items())[:500])

    return TaxonResult(results=TaxonDict(taxon=results), description=desc, count = len(results))

app = VersionedFastAPI(app, enable_latest=True)
