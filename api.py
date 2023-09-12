from typing import Dict, List

from fastapi import FastAPI
from fastapi_versioning import VersionedFastAPI, version
from pydantic import BaseModel

from model import DataModel

description = """
LOTUSFast API helps you do awesome stuff. 🚀
"""


class Item(BaseModel):
    structure_wid: int | None = None  # 3613679
    molecule: str | None = None  # "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O"
    substructure_search: bool | None = None  # False
    similarity_level: float | None = None  # 0.8
    taxon_wid: int | None = None  # 158572
    taxon_name: str | None = None  # "Gentiana lutea"

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

class CoupleDict(BaseModel):
    couple: Dict

class CoupleResult(BaseModel):
    results: CoupleDict
    count: int
    description: str

class MainResult(BaseModel):
    results_couples: CoupleDict
    results_structures: StructureDict
    results_taxa: TaxonDict


def create_main(
        dm: DataModel,
        structure_wid: int = None,
        molecule: str = None,
        substructure_search: bool = None,
        similarity_level: float = None,
        taxon_wid: int = None,
        taxon_name: str = None
) -> MainResult:
    
    ## Structural part
    structures = dm.get_compounds()
    if molecule:
        ids = list(dm.compound_search(molecule)) if not substructure_search else list(
            dm.compound_search_substructure(molecule))
        structures = [id for id, score in ids if score == 1] if not substructure_search and not similarity_level else [
            id for id, score in ids if
            score >= similarity_level]
    if structure_wid:
        structures = [x for x in structures if x == structure_wid]
    tax_from_str = [item for wid in structures for item in dm.get_taxa_containing_compound(wid)]

    ## Taxal part
    taxa = dm.get_taxa()
    ## TODO ADD
    # get_taxonomic_tree
    if taxon_name:
        ## Not using resolve_taxon for now
        ids = list(dm.get_taxa_with_name_containing(taxon_name))
        taxon_names = dm.get_taxon_name_from_list_of_wid(ids)
        taxa = dict(zip(ids, taxon_names))
    if taxon_wid:
        taxa = {key: value for key, value in taxa.items() if key == taxon_wid}
    str_from_tax = [item for wid in taxa for item in dm.get_compounds_of_taxon(wid)]

    ## Filter both ways
    structures = [value for value in structures if value in str_from_tax]
    taxa = {key: value for key, value in taxa.items() if key in tax_from_str}

    ## For dev
    structures = structures[:500]
    smiles = dm.get_compound_smiles_from_list_of_wid(structures)
    taxa = dict(list(taxa.items())[:500])

    ## Couples part
    couples = {wid: [] for wid in taxa for item in dm.get_compounds_of_taxon(wid)}
    for wid, item in [(wid, item) for wid in taxa for item in dm.get_compounds_of_taxon(wid)]:
        couples[wid].append(item)

    ## Final
    s = StructureDict(structure_id=structures, structure_smiles=smiles)
    t = TaxonDict(taxon=taxa)
    c = CoupleDict(couple=couples)

    return MainResult(results_couples=c, results_structures=s, results_taxa=t)

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
async def create_couples(item: Item) -> CoupleResult:
    results = create_main(
        dm=dm,
        structure_wid=item.structure_wid,
        molecule=item.molecule,
        substructure_search=item.substructure_search,
        similarity_level=item.similarity_level,
        taxon_wid=item.taxon_wid,
        taxon_name=item.taxon_name)

    return CoupleResult(results=results.results_couples,
                        description="Couples matching the query",
                        count=sum(len(items) for items in results.results_couples.couple.values()))

@app.post("/structures/")
@version(1, 0)
async def create_structures(item: Item) -> StructureResult:
    results = create_main(
        dm=dm,
        structure_wid=item.structure_wid,
        molecule=item.molecule,
        substructure_search=item.substructure_search,
        similarity_level=item.similarity_level,
        taxon_wid=item.taxon_wid,
        taxon_name=item.taxon_name)

    return StructureResult(results=results.results_structures,
                           description="Structures matching the query",
                           count=len(results.results_structures.structure_id))

@app.post("/taxa/")
@version(1, 0)
async def create_taxa(item: Item) -> TaxonResult:
    results = create_main(
        dm=dm,
        structure_wid=item.structure_wid,
        molecule=item.molecule,
        substructure_search=item.substructure_search,
        similarity_level=item.similarity_level,
        taxon_wid=item.taxon_wid,
        taxon_name=item.taxon_name)

    return TaxonResult(results=results.results_taxa,
                       description="Taxa matching the query",
                       count=len(results.results_taxa.taxon.items()))

app = VersionedFastAPI(app, enable_latest=True)
