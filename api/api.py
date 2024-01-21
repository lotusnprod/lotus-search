import logging

from fastapi import Depends, FastAPI
from fastapi_versioning import VersionedFastAPI, version

from api.models import (
    Item,
    ReferenceInfo,
    ReferenceResult,
    StructureInfo,
    StructureResult,
    TaxonInfo,
    TaxonResult,
    TripletResult,
)
from api.queries import (
    get_references_for_item,
    get_structures_for_item,
    get_taxa_for_item,
    get_triplets_for_item,
)
from model.data_model import DataModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

description = """
LOTUSFast API helps you do awesome stuff. 🚀
"""


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
async def search_triplets(
    item: Item, dm: DataModel = Depends(DataModel)
) -> TripletResult:
    triplets = get_triplets_for_item(item, dm)

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
async def search_structures(
    item: Item, dm: DataModel = Depends(DataModel)
) -> StructureResult:
    dict_items = get_structures_for_item(item, dm)

    return StructureResult(
        ids=dict_items.keys(),
        structures={
            sid: StructureInfo(smiles=value) for sid, value in dict_items.items()
        },
        description="Structures matching the query",
        count=len(dict_items),
    )


@app.post("/taxa")
@version(1, 0)
async def search_taxa(item: Item, dm: DataModel = Depends(DataModel)) -> TaxonResult:
    dict_items = get_taxa_for_item(item, dm)

    return TaxonResult(
        ids=dict_items.keys(),
        taxa={tid: TaxonInfo(name=value) for tid, value in dict_items.items()},
        description="Taxa matching the query",
        count=len(dict_items),
    )


@app.post("/references")
@version(1, 0)
async def search_references(
    item: Item, dm: DataModel = Depends(DataModel)
) -> ReferenceResult:
    dict_items = get_references_for_item(item, dm)

    return ReferenceResult(
        ids=dict_items.keys(),
        references={rid: ReferenceInfo(doi=value) for rid, value in dict_items.items()},
        description="References matching the query",
        count=len(dict_items),
    )


app = VersionedFastAPI(app, enable_latest=True)
