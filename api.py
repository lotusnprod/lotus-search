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
    structure_wid: Annotated[
        int | None,
        Query(
            alias="structure_wid",
            description="Wikidata identifier of the structure (without the Q).",
            example="3613679",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    molecule: Annotated[
        str | None,
        Query(
            alias="molecule",
            description="A MOL file or SMILES of the structure of the structure or part of it.",
            example="C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
        ),
    ] = None,
        substructure_search: Annotated[
            bool | None,
            Query(
                alias="substructure_search",
                description="Search by substructure.",
                example="false"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    similarity_level: Annotated[
        float | None,
        Query(
            alias="similarity_level",
            description="Similarity level cut-off (basic tanimoto-like search). Does nothing is substructure_search is true.",
            example="0.8"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    taxon_wid: Annotated[
        int | None,
        Query(
            alias="taxon_wid",
            description="Wikidata identifier of the taxon (without the Q).",
            example="158572",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    taxon_name: Annotated[
        str | None,
        Query(
            alias="taxon_name",
            description="The name searched (can be partial and slightly incorrect).",
            example="Gentiana lutea"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None
) -> CoupleResult:
    desc = "Couples matching the query"
    ## TODO
    return CoupleResult(results=CoupleDict(structure_id=couples, structure_smiles=smiles, taxon=taxa), description=desc, count = len(couples))

@app.post("/structures/")
@version(1, 0)
async def create_structures(
    structure_wid: Annotated[
        int | None,
        Query(
            alias="structure_wid",
            description="Wikidata identifier of the structure (without the Q).",
            example="3613679",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    molecule: Annotated[
        str | None,
        Query(
            alias="molecule",
            description="A MOL file or SMILES of the structure of the structure or part of it.",
            example="C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
        ),
    ] = None,
        substructure_search: Annotated[
            bool | None,
            Query(
                alias="substructure_search",
                description="Search by substructure.",
                example="false"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    similarity_level: Annotated[
        float | None,
        Query(
            alias="similarity_level",
            description="Similarity level cut-off (basic tanimoto-like search). Does nothing is substructure_search is true.",
            example="0.8"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    taxon_wid: Annotated[
        int | None,
        Query(
            alias="taxon_wid",
            description="Wikidata identifier of the taxon (without the Q).",
            example="158572",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    taxon_name: Annotated[
        str | None,
        Query(
            alias="taxon_name",
            description="The name searched (can be partial and slightly incorrect).",
            example="Gentiana lutea"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None
) -> StructureResult:
    desc = "Structures matching the query"
    ## TODO add taxal filtering
    structures = dm.get_compounds()
    if molecule:
        ids =  list(dm.compound_search(molecule)) if not substructure_search else list(dm.compound_search_substructure(molecule))
        structures = [id for id, score in ids if score == 1] if not substructure_search and not similarity_level else [id for id, score in ids if
                                                                                        score >= similarity_level]
    if structure_wid:
        structures = [x for x in structures if x == structure_wid]
    ## For dev
    structures = structures[:500]
    smiles = dm.get_compound_smiles_from_list_of_wid(structures)

    return StructureResult(results=StructureDict(structure_id=structures, structure_smiles=smiles), description=desc, count = len(structures))

@app.post("/taxa/")
@version(1, 0)
async def create_taxa(
        structure_wid: Annotated[
        int | None,
        Query(
            alias="structure_wid",
            description="Wikidata identifier of the structure (without the Q).",
            example="3613679",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    molecule: Annotated[
        str | None,
        Query(
            alias="molecule",
            description="A MOL file or SMILES of the structure of the structure or part of it.",
            example="C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
        ),
    ] = None,
        substructure_search: Annotated[
            bool | None,
            Query(
                alias="substructure_search",
                description="Search by substructure.",
                example="false"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    similarity_level: Annotated[
        float | None,
        Query(
            alias="similarity_level",
            description="Similarity level cut-off (basic tanimoto-like search). Does nothing is substructure_search is true.",
            example="0.8"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    taxon_wid: Annotated[
        int | None,
        Query(
            alias="taxon_wid",
            description="Wikidata identifier of the taxon (without the Q).",
            example="158572",
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None,
    taxon_name: Annotated[
        str | None,
        Query(
            alias="taxon_name",
            description="The name searched (can be partial and slightly incorrect).",
            example="Gentiana lutea"
            # min_length=3,
            # max_length=50,
            # pattern="^fixedquery$"
            ),
    ] = None
) -> TaxonResult:
    desc = "Taxa matching the query"
    ## TODO add structural filtering
    taxa = dm.get_taxa()
    if taxon_name:
        ids = list(dm.get_taxa_with_name_containing(taxon_name))
        taxon_names = dm.get_taxon_name_from_list_of_wid(ids)
        taxa = dict(zip(ids, taxon_names))

    if taxon_wid:
        taxa = {key: value for key, value in taxa.items() if key == taxon_wid}
    ## For dev
    taxa = dict(list(taxa.items())[:500])

    return TaxonResult(results=TaxonDict(taxon=taxa), description=desc, count = len(taxa))

app = VersionedFastAPI(app, enable_latest=True)
