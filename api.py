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
    similarity_level: float = 1.0  # exact if not given
    taxon_wid: int | None = None  # 158572
    taxon_name: str | None = None  # "Gentiana lutea"


class StructureResult(BaseModel):
    results: dict[int, str]  ## AR: Later you change this str by a StructureInfo (just make the method in the model)
    count: int
    description: str


class TaxonResult(BaseModel):
    results: dict[int, str]  ## AR: Later you change this str by a TaxonInfo (just make the method in the model)
    count: int
    description: str


class CoupleResult(BaseModel):
    results: list[dict[str, int]]
    count: int
    description: str


dm = DataModel()
## FOR AR: We use that all the time lets precalculate. Making it a set for the next steps
## Should likely move in the data model if that's used all the time
all_structures = set(dm.get_compounds())
all_taxa = dm.get_taxa()


def get_matching_structures_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns None if the item do not filter by structures"""
    if item.molecule is None and item.structure_wid is None:
        return None

    ## FOR AR: I would argue here that if the user gives a WID, we shouldn't do the structure search at all
    ## And this needs to be explained in the API doc
    if item.structure_wid:
        if item.structure_wid in all_structures:
            return {item.structure_wid}
        else:
            return {}

    if item.molecule:
        if item.substructure_search:  ## FOR AR: Make this kind of IFs on multiple lines as it's easier to read and understand
            results = dm.compound_search_substructure(item.molecule)  ## FOR AR: This was already a list
            structures = {_id for _id, _ in results}
        else:
            results = dm.compound_search(item.molecule)  ## FOR AR: Same
            structures = {_id for _id, score in results if score >= item.similarity_level}  ## FOR AR: ID is a reserved keyword in python
    else:
        structures = all_structures  ## FOR AR: Load them all only if you need them all, and we don't modify that object content ever

    return structures

def get_taxa_matching_query(dm: DataModel, item: Item) -> set[int] | None:
    """If it returns None, then no element is filtering by taxa, otherwise it returns the WID of matching taxa"""
    if not item.taxon_wid and not item.taxon_name:
        return None

    taxa = set()

    if item.taxon_wid:  ## AR: If we give a taxon_wid it takes priority, that's probably the right thing to do
        taxa = {item.taxon_wid}
    else:
        if item.taxon_name:
            taxa = set(dm.get_taxa_with_name_containing(item.taxon_name))

    return taxa

def get_matching_structures_from_taxa_in_item(dm: DataModel, item: Item) -> set[int] | None:
    # We need to get all the matching taxa
    taxa = get_taxa_matching_query(dm, item)

    if taxa is None:  # We are not filtering by taxon at all
        return None

    out = {dm.get_compounds_of_taxon(taxon) for taxon in taxa} # AR: We could have a parameter "recursive" in the query to have all the compounds from the parents too

    return out

def get_matching_taxa_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    # We need to get all the matching structures
    structures = get_matching_structures_from_structure_in_item(dm, item)

    if structures is None:  # We are not filtering by structure at all
        return None

    out = {dm.get_taxa_containing_compound(structure) for structure in structures}

    return out


app = FastAPI(title="LOTUS FastAPI", description=description,
              summary="An awesome way to access natural products related data.")


@app.post("/couples")
@version(1, 0)
async def create_couples(item: Item) -> CoupleResult:
    selected_structures = get_matching_structures_from_structure_in_item(dm, item)
    selected_taxa = get_taxa_matching_query(dm, item)
    couples = set()

    if selected_taxa is not None:
        for taxon in selected_taxa:
            structures = dm.get_compounds_of_taxon(taxon)  # Same thing could add recursive here, and we should rename this to get_structures_of_taxon
            # To add reference, just add it to the if here
            couples.update({(taxon, structure) for structure in structures if structure in selected_structures})

    ## AR: Left as an exercise to add taxon name or molecule details or what not
    ## Maybe to not resolve by default, and let them do the complex query if they need to? Not everybody need names and
    ## that make everything much slower
    return CoupleResult(results=[{"taxon_wid": taxon, "structure_wid": structure} for taxon, structure in couples],
                        description="Couples matching the query", count=len(couples))


@app.post("/structures")
@version(1, 0)
async def create_structures(item: Item) -> StructureResult:
    # We want the set of all the structures that match a query, or the specific structure given
    selected_structures = get_matching_structures_from_structure_in_item(dm, item)

    # We want the set of all the molecules found in the given taxa
    matching_structures_by_taxa = get_matching_structures_from_taxa_in_item(dm, item)

    # We want the intersection of both (and we can do the same for the references later)

    ## AR: The other solution is to return all taxa or all compounds when the query doesn't match anything
    ## This should have no impact on performance and may make that easier to read
    if matching_structures_by_taxa is not None:
        if selected_structures is None:
            selected_structures = matching_structures_by_taxa
        else:
            selected_structures = selected_structures & matching_structures_by_taxa

    if selected_structures is None:
        selected_structures = set()

    return StructureResult(results=dm.get_dict_of_wid_to_smiles(selected_structures),
                           description="Structures matching the query",
                           count=len(selected_structures))


@app.post("/taxa")
@version(1, 0)
async def create_taxa(item: Item) -> TaxonResult:
    # We want the set of all the taxa matching the query
    selected_taxa = get_taxa_matching_query(dm, item)

    # We want the set of all the taxa which have molecules matching the query
    matching_taxa_by_structure = get_matching_taxa_from_structure_in_item(dm, item)
    # We want the intersection of both (and we can do the same for the references later)
    if matching_taxa_by_structure is not None:
        if selected_taxa is None:
            selected_taxa = matching_taxa_by_structure
        else:
            selected_taxa = selected_taxa & matching_taxa_by_structure

    if selected_taxa is None:
        selected_taxa = set()

    return TaxonResult(results=dm.get_dict_of_wid_to_taxon_name(selected_taxa),
                       description="Taxa matching the query",
                       count=len(selected_taxa))

app = VersionedFastAPI(app, enable_latest=True)
