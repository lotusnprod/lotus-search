import logging

from fastapi import FastAPI, HTTPException
from fastapi_versioning import VersionedFastAPI, version

from api.models import (CoupleResult, Item, StructureInfo, StructureResult,
                        TaxonInfo, TaxonResult)
from model import DataModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


description = """
LOTUSFast API helps you do awesome stuff. 🚀
"""

dm = DataModel()
# Should likely move in the data model if that's used all the time

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
    },
)


def get_matching_structures_from_structure_in_item(
    dm: DataModel, item: Item
) -> set[int]:
    """Returns all_structures if the item do not filter by structure, else returns the WID of matching structures"""
    if item.structure is None and item.structure_wid is None:
        return dm.structures_set()
    elif item.structure and item.structure_wid:
        raise HTTPException(
            status_code=500,
            detail=f"You cannot provide both 'structure' and 'structure_wid'",
        )
    else:
        # This needs to be explained in the API doc
        if item.structure_wid:
            if item.structure_wid in dm.structures_set():
                return {item.structure_wid}
        else:
            if item.structure:
                if item.substructure_search:
                    try:
                        results = dm.structure_search_substructure(item.structure)
                        structures = {_id for _id, _ in results}
                    except ValueError:
                        raise HTTPException(
                            status_code=500,
                            detail=f"The structure given is invalid: {item.structure}",
                        )
                else:
                    try:
                        results = dm.structure_search(item.structure)
                        structures = {
                            _id
                            for _id, score in results
                            if score >= item.similarity_level
                        }
                    except ValueError:
                        raise HTTPException(
                            status_code=500,
                            detail=f"The structure given is invalid: {item.structure}",
                        )
            else:
                structures = dm.structures_set()

            return structures


def get_matching_taxa_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns all_taxa if the item do not filter by taxon, else returns the WID of matching taxa or None if no taxa requested"""
    if item.taxon_wid is None and item.taxon_name is None:
        return None
    else:
        # This needs to be explained in the API doc
        if item.taxon_wid:
            if item.taxon_wid in dm.get_taxa():
                return {item.taxon_wid}
        else:
            if item.taxon_name:
                taxa = set(dm.get_taxa_with_name_containing(item.taxon_name))
                if taxa is None:
                    taxa = dm.get_taxa()
            else:
                taxa = dm.get_taxa()

        return taxa


def get_matching_structures_from_taxon_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching taxa
    taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return None

    # We could have a parameter "recursive" in the query to have all the structures from the parents too
    out = set()
    for taxon in taxa:
        out.update(dm.get_structures_of_taxon(taxon))

    return out


def get_matching_taxa_from_structure_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching structures
    structures = get_matching_structures_from_structure_in_item(dm, item)

    out = set()
    for structure in structures:
        out.update(dm.get_taxa_containing_structure(structure))

    return out


@app.post("/couples")
@version(1, 0)
async def search_couples(item: Item) -> CoupleResult:
    selected_structures = get_matching_structures_from_structure_in_item(dm, item)
    selected_taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    structures_of_selected_taxa = {
        taxon: dm.get_structures_of_taxon(taxon) for taxon in selected_taxa
    }

    couples = {
        (structure, taxon)
        for taxon, structures in structures_of_selected_taxa.items()
        for structure in structures
        if structure in selected_structures
    }

    if item.limit == 0:
        couples = list(couples)
    else:
        couples = list(couples)[: item.limit]

    return CoupleResult(
        ids=[{"structure": structure, "taxon": taxon} for structure, taxon in couples],
        structures={
            wid: StructureInfo(smiles=value)
            for wid, value in dm.get_dict_of_wid_to_smiles(
                [first_value for first_value, _ in couples]
            ).items()
        },
        taxa={
            wid: TaxonInfo(name=value)
            for wid, value in dm.get_dict_of_wid_to_taxon_name(
                [taxon_name for _, taxon_name in couples]
            ).items()
        },
        description="Couples matching the query",
        count=len(couples),
    )


@app.post("/structures")
@version(1, 0)
async def search_structures(item: Item) -> StructureResult:
    # We want the set of all the structures that match a query, or the specific structure given
    matching_structures_by_structure = get_matching_structures_from_structure_in_item(
        dm, item
    )
    # We want the set of all the structures found in the given taxa
    matching_structures_by_taxon = get_matching_structures_from_taxon_in_item(dm, item)

    # We want the intersection of both (and we can do the same for the references later)
    # But if one of the sets is fully empty
    if matching_structures_by_taxon is None:
        matching_structures = matching_structures_by_structure
    else:
        matching_structures = (
            matching_structures_by_structure & matching_structures_by_taxon
        )

    if item.limit == 0:
        items = dm.get_dict_of_wid_to_smiles(matching_structures).items()
    else:
        items = dm.get_dict_of_wid_to_smiles(matching_structures).items()[: item.limit]

    return StructureResult(
        ids=matching_structures,
        structures={wid: StructureInfo(smiles=value) for wid, value in items},
        description="Structures matching the query",
        count=len(matching_structures),
    )


@app.post("/taxa")
@version(1, 0)
async def search_taxa(item: Item) -> TaxonResult:
    # We want the set of all the taxa matching the query
    matching_taxa_by_taxon = get_matching_taxa_from_taxon_in_item(dm, item)

    # We want the set of all the taxa which have structures matching the query
    matching_taxa_by_structure = get_matching_taxa_from_structure_in_item(dm, item)

    # We want the intersection of both (and we can do the same for the references later)
    matching_taxa = matching_taxa_by_taxon & matching_taxa_by_structure

    if item.limit == 0:
        items = dm.get_dict_of_wid_to_taxon_name(matching_taxa).items()
    else:
        items = dm.get_dict_of_wid_to_taxon_name(matching_taxa).items()[: item.limit]

    return TaxonResult(
        ids=matching_taxa,
        taxa={wid: TaxonInfo(name=value) for wid, value in items},
        description="Taxa matching the query",
        count=len(matching_taxa),
    )


app = VersionedFastAPI(app, enable_latest=True)
