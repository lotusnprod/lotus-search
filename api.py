from fastapi import FastAPI, HTTPException
from fastapi_versioning import VersionedFastAPI, version

from model import CoupleResult, DataModel, Item, StructureInfo, StructureResult, TaxonInfo, TaxonResult

description = """
LOTUSFast API helps you do awesome stuff. 🚀
"""

dm = DataModel()
# Should likely move in the data model if that's used all the time
all_structures = set(dm.get_compounds())
all_taxa = dm.get_taxa()

app = FastAPI(
    title="LOTUS FastAPI",
    description=description,
    summary="An awesome way to access natural products related data.",
    version="1.0",
    # TODO
    # terms_of_service="http://example.com/terms/",
    # contact={
    #     "name": "Deadpoolio the Amazing",
    #     "url": "http://x-force.example.com/contact/",
    #     "email": "dp@x-force.example.com",
    # },
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    }
)

def get_matching_structures_from_structure_in_item(dm: DataModel, item: Item) -> set[int]:
    """Returns all_structures if the item do not filter by structure, else returns the WID of matching structures"""
    if item.molecule is None and item.structure_wid is None:
        return all_structures
    else:
        # This needs to be explained in the API doc
        if item.structure_wid:
            if item.structure_wid in all_structures:
                return {item.structure_wid}
        else:
            if item.molecule:
                if item.substructure_search:
                    try:
                        results = dm.compound_search_substructure(item.molecule)
                        structures = {_id for _id, _ in results}
                    except ValueError:
                        raise HTTPException(status_code=500, detail=f"The structure given is invalid: {item.molecule}")
                else:
                    try:
                        results = dm.compound_search(item.molecule)
                        structures = {_id for _id, score in results if score >= item.similarity_level}
                    except ValueError:
                        raise HTTPException(status_code=500, detail=f"The structure given is invalid: {item.molecule}")
            else:
                structures = all_structures

            return structures


def get_matching_taxa_from_taxon_in_item(dm: DataModel, item: Item) -> set[int]:
    """Returns all_taxa if the item do not filter by taxon, else returns the WID of matching taxa"""
    if item.taxon_wid is None and item.taxon_name is None:
        return all_taxa
    else:
        # This needs to be explained in the API doc
        if item.taxon_wid:
            if item.taxon_wid in all_taxa:
                return {item.taxon_wid}
        else:
            if item.taxon_name:
                taxa = set(dm.get_taxa_with_name_containing(item.taxon_name))
                if taxa is None:
                    taxa = all_taxa
            else:
                taxa = all_taxa

        return taxa


def get_matching_structures_from_taxon_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching taxa
    taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    # We could have a parameter "recursive" in the query to have all the compounds from the parents too
    out = set()
    for taxon in taxa:
        out.update(dm.get_compounds_of_taxon(taxon))

    return out


def get_matching_taxa_from_structure_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching structures
    structures = get_matching_structures_from_structure_in_item(dm, item)

    out = set()
    for structure in structures:
        out.update(dm.get_taxa_containing_compound(structure))

    return out


@app.post("/couples")
@version(1, 0)
async def search_couples(item: Item) -> CoupleResult:
    selected_structures = get_matching_structures_from_structure_in_item(dm, item)
    selected_taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    compounds_of_selected_taxa = {
        taxon: dm.get_compounds_of_taxon(taxon) for taxon in selected_taxa
    }

    couples = {
        (structure, taxon)
        for taxon, structures in compounds_of_selected_taxa.items()
        for structure in structures
        if structure in selected_structures
    }

    return CoupleResult(
        ids=[{"structure": structure, "taxon": taxon} for structure, taxon in couples],
        structures={wid: StructureInfo(smiles=value) for wid, value in
                    dm.get_dict_of_wid_to_smiles([first_value for first_value, _ in couples]).items()},
        taxa={wid: TaxonInfo(name=value) for wid, value in
              dm.get_dict_of_wid_to_taxon_name([taxon_name for _, taxon_name in couples]).items()},
        description="Couples matching the query",
        count=len(couples))


@app.post("/structures")
@version(1, 0)
async def search_structures(item: Item) -> StructureResult:
    # We want the set of all the structures that match a query, or the specific structure given
    matching_structures_by_structure = get_matching_structures_from_structure_in_item(dm, item)
    # We want the set of all the molecules found in the given taxa
    matching_structures_by_taxon = get_matching_structures_from_taxon_in_item(dm, item)

    # We want the intersection of both (and we can do the same for the references later)
    matching_structures = matching_structures_by_structure & matching_structures_by_taxon
    return StructureResult(
        ids=matching_structures,
        structures={wid: StructureInfo(smiles=value) for wid, value in
                    dm.get_dict_of_wid_to_smiles(matching_structures).items()},
        description="Structures matching the query",
        count=len(matching_structures)
    )


@app.post("/taxa")
@version(1, 0)
async def search_taxa(item: Item) -> TaxonResult:
    # We want the set of all the taxa matching the query
    matching_taxa_by_taxon = get_matching_taxa_from_taxon_in_item(dm, item)

    # We want the set of all the taxa which have molecules matching the query
    matching_taxa_by_structure = get_matching_taxa_from_structure_in_item(dm, item)

    # We want the intersection of both (and we can do the same for the references later)
    matching_taxa = matching_taxa_by_taxon & matching_taxa_by_structure

    return TaxonResult(
        ids=matching_taxa,
        infos={wid: TaxonInfo(name=value) for wid, value in dm.get_dict_of_wid_to_taxon_name(matching_taxa).items()},
        description="Taxa matching the query",
        count=len(matching_taxa)
    )

app = VersionedFastAPI(app, enable_latest=True)
