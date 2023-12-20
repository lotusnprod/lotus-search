import logging

from fastapi import FastAPI, HTTPException
from fastapi_versioning import VersionedFastAPI, version

from api.models import (CoupleResult, Item, ReferenceInfo, ReferenceResult,
                        StructureInfo, StructureResult, TaxonInfo, TaxonResult)
from api.queries import (  # get_matching_references_from_couple_in_item,
    get_matching_references_from_reference_in_item,
    get_matching_references_from_structure_in_item,
    get_matching_references_from_taxon_in_item,
    get_matching_structures_from_reference_in_item,
    get_matching_structures_from_structure_in_item,
    get_matching_structures_from_taxon_in_item,
    get_matching_taxa_from_reference_in_item,
    get_matching_taxa_from_structure_in_item,
    get_matching_taxa_from_taxon_in_item)
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


@app.post("/couples")
@version(1, 0)
async def search_couples(item: Item) -> CoupleResult:
    selected_structures = get_matching_structures_from_structure_in_item(dm, item)
    selected_taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    structures_of_selected_taxa = (
        {taxon: dm.get_structures_of_taxon(taxon) for taxon in selected_taxa}
        if selected_taxa is not None
        else {
            taxon: dm.get_structures_of_taxon(taxon)
            for structure in selected_structures
            for taxon in dm.get_taxa_containing_structure(structure)
        }
    )

    # TODO add references

    couples = {
        (structure, taxon)
        for taxon, structures in structures_of_selected_taxa.items()
        if selected_structures is not None
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
            for wid, value in dm.get_dict_of_sid_to_smiles(
                [first_value for first_value, _ in couples]
            ).items()
        },
        taxa={
            wid: TaxonInfo(name=value)
            for wid, value in dm.get_dict_of_tid_to_taxon_name(
                [taxon_name for _, taxon_name in couples]
            ).items()
        },
        description="Couples matching the query",
        count=len(couples),
    )


@app.post("/structures")
@version(1, 0)
async def search_structures(item: Item) -> StructureResult:
    matching_structures_by_structure = get_matching_structures_from_structure_in_item(
        dm, item
    )
    matching_structures_by_taxon = get_matching_structures_from_taxon_in_item(dm, item)
    matching_structures_by_reference = get_matching_structures_from_reference_in_item(
        dm, item
    )

    non_empty_sets = [
        s
        for s in [
            matching_structures_by_reference,
            matching_structures_by_taxon,
            matching_structures_by_structure,
        ]
        if s
    ]
    matching_structures = set.intersection(*non_empty_sets) if non_empty_sets else set()

    items = list(dm.get_dict_of_sid_to_smiles(matching_structures).items())

    if item.limit == 0:
        items = items
    else:
        items = items[: item.limit]

    return StructureResult(
        ids=matching_structures,
        structures={sid: StructureInfo(smiles=value) for sid, value in items},
        description="Structures matching the query",
        count=len(matching_structures),
    )


@app.post("/taxa")
@version(1, 0)
async def search_taxa(item: Item) -> TaxonResult:
    matching_taxa_by_taxon = get_matching_taxa_from_taxon_in_item(dm, item)
    matching_taxa_by_structure = get_matching_taxa_from_structure_in_item(dm, item)
    matching_taxa_by_reference = get_matching_taxa_from_reference_in_item(dm, item)

    non_empty_sets = [
        s
        for s in [
            matching_taxa_by_reference,
            matching_taxa_by_structure,
            matching_taxa_by_taxon,
        ]
        if s
    ]
    matching_taxa = set.intersection(*non_empty_sets) if non_empty_sets else set()

    items = list(dm.get_dict_of_tid_to_taxon_name(matching_taxa).items())

    if item.limit == 0:
        items = items
    else:
        items = items[: item.limit]

    return TaxonResult(
        ids=matching_taxa,
        taxa={tid: TaxonInfo(name=value) for tid, value in items},
        description="Taxa matching the query",
        count=len(matching_taxa),
    )


@app.post("/references")
@version(1, 0)
async def search_references(item: Item) -> ReferenceResult:
    matching_references_by_reference = get_matching_references_from_reference_in_item(
        dm, item
    )
    # matching_references_by_couple = get_matching_references_from_couple_in_item(
    #     dm, item
    # )
    matching_references_by_structure = get_matching_references_from_structure_in_item(
        dm, item
    )
    matching_references_by_taxon = get_matching_references_from_taxon_in_item(dm, item)

    non_empty_sets = [
        s
        for s in [
            matching_references_by_reference,
            matching_references_by_structure,
            matching_references_by_taxon,
        ]
        if s
    ]
    matching_references = set.intersection(*non_empty_sets) if non_empty_sets else set()

    items = list(dm.get_dict_of_rid_to_reference_doi(matching_references).items())

    if item.limit == 0:
        items = items
    else:
        items = items[: item.limit]

    return ReferenceResult(
        ids=matching_references,
        references={rid: ReferenceInfo(doi=value) for rid, value in items},
        description="References matching the query",
        count=len(matching_references),
    )


app = VersionedFastAPI(app, enable_latest=True)
