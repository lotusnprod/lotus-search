import logging

from fastapi import Depends, FastAPI
from fastapi_versioning import VersionedFastAPI, version

from api.models import (Item, ReferenceInfo, ReferenceResult,
                        StructureInfo, StructureResult, TaxonInfo, TaxonResult, TripletResult)
from api.queries import (combine_and_filter_outputs, get_matching_references_from_reference_in_item,
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


def get_dm():
    return DataModel()


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


@app.post("/triplets")
@version(1, 0)
async def search_triplets(item: Item, dm: DataModel = Depends(get_dm)) -> TripletResult:
    selected_references = get_matching_references_from_reference_in_item(dm, item)
    selected_structures = get_matching_structures_from_structure_in_item(dm, item)
    selected_taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    triplets_set = dm.get_triples_for(reference_ids=selected_references,
                                      structure_ids=selected_structures,
                                      taxon_ids=selected_taxa)

    if item.limit == 0:
        triplets = list(triplets_set)
    else:
        triplets = list(triplets_set)[:item.limit]

    return TripletResult(
        triplets=triplets,
        references={
            wid: ReferenceInfo(doi=value)
            for wid, value in dm.get_dict_of_rid_to_reference_doi(
                [reference_id for reference_id, _, _ in triplets]
            ).items()
        },
        structures={
            wid: StructureInfo(smiles=value)
            for wid, value in dm.get_dict_of_sid_to_smiles(
                [structure_id for _, structure_id, _ in triplets]
            ).items()
        },
        taxa={
            wid: TaxonInfo(name=value)
            for wid, value in dm.get_dict_of_tid_to_taxon_name(
                [taxon_id for _, _, taxon_id in triplets]
            ).items()
        },
        description="Triplets matching the query",
        count=len(triplets),
    )


@app.post("/structures")
@version(1, 0)
async def search_structures(item: Item, dm: DataModel = Depends(get_dm)) -> StructureResult:
    matching_structures_by_structure = get_matching_structures_from_structure_in_item(
        dm, item
    )
    matching_structures_by_taxon = get_matching_structures_from_taxon_in_item(dm, item)
    matching_structures_by_reference = get_matching_structures_from_reference_in_item(
        dm, item
    )

    ids = combine_and_filter_outputs([matching_structures_by_reference,
                                      matching_structures_by_taxon,
                                      matching_structures_by_structure], limit=item.limit)

    dict_items = dm.get_dict_of_sid_to_smiles(ids)

    return StructureResult(
        ids=dict_items.keys(),
        structures={sid: StructureInfo(smiles=value) for sid, value in dict_items.items()},
        description="Structures matching the query",
        count=len(dict_items),
    )


@app.post("/taxa")
@version(1, 0)
async def search_taxa(item: Item, dm: DataModel = Depends(get_dm)) -> TaxonResult:
    matching_taxa_by_taxon = get_matching_taxa_from_taxon_in_item(dm, item)
    matching_taxa_by_structure = get_matching_taxa_from_structure_in_item(dm, item)
    matching_taxa_by_reference = get_matching_taxa_from_reference_in_item(dm, item)

    ids = combine_and_filter_outputs([matching_taxa_by_reference,
                                      matching_taxa_by_structure,
                                      matching_taxa_by_taxon], limit=item.limit)

    dict_items = dm.get_dict_of_tid_to_taxon_name(ids)

    return TaxonResult(
        ids=dict_items.keys(),
        taxa={tid: TaxonInfo(name=value) for tid, value in dict_items.items()},
        description="Taxa matching the query",
        count=len(dict_items),
    )


@app.post("/references")
@version(1, 0)
async def search_references(item: Item, dm: DataModel = Depends(get_dm)) -> ReferenceResult:
    matching_references_by_reference = get_matching_references_from_reference_in_item(
        dm, item
    )
    matching_references_by_structure = get_matching_references_from_structure_in_item(
        dm, item
    )
    matching_references_by_taxon = get_matching_references_from_taxon_in_item(dm, item)

    ids = combine_and_filter_outputs([matching_references_by_reference,
                                      matching_references_by_structure,
                                      matching_references_by_taxon], limit=item.limit)

    dict_items = dm.get_dict_of_rid_to_reference_doi(ids)

    return ReferenceResult(
        ids=dict_items.keys(),
        references={rid: ReferenceInfo(doi=value) for rid, value in dict_items.items()},
        description="References matching the query",
        count=len(dict_items),
    )


app = VersionedFastAPI(app, enable_latest=True)
