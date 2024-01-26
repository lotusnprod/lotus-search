import logging
from typing import Any

from fastapi import HTTPException

from api.models import (
    FilterItem,
    Item,
    ReferenceItem,
    StructureFilterItem,
    StructureItem,
    TaxonFilterItem,
    TaxonItem,
)
from model.data_model import DataModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def references_from_reference_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns the WID of matching references."""
    wid = item.reference.wid
    doi = item.reference.doi

    if doi and wid:
        raise HTTPException(
            status_code=500,
            detail=f"You cannot provide both 'doi' and 'wid'",
        )
    elif wid is not None or doi is not None:
        # This needs to be explained in the API doc
        if wid:
            return dm.get_reference_with_id(wid)
        return dm.get_references_with_doi(doi)

    return None


def structures_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns the WID of matching structures."""
    structures = None

    wid = item.structure.wid
    molecule = item.structure.molecule
    sub = item.filter.structure.substructure_search
    sim = item.filter.structure.similarity_level

    if molecule and wid:
        raise HTTPException(
            status_code=500,
            detail=f"You cannot provide both 'molecule' and 'wid'",
        )
    elif molecule is not None or wid is not None:
        # This needs to be explained in the API doc
        if wid:
            if wid in dm.structures_set():
                return {wid}
            else:
                return set()

        if sub:
            try:
                results = dm.structure_search_substructure(molecule)
                structures = {_id for _id, _ in results}
            except ValueError:
                raise HTTPException(
                    status_code=500,
                    detail=f"The structure given is invalid: {molecule}",
                )
        else:
            try:
                results = dm.structure_search(molecule)
                structures = {_id for _id, score in results if score >= sim}
            except ValueError:
                raise HTTPException(
                    status_code=500,
                    detail=f"The structure given is invalid: {molecule}",
                )
        return structures

    return structures


def taxa_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns the WID of matching taxa."""
    wid = item.taxon.wid
    name = item.taxon.name

    if name and wid:
        raise HTTPException(
            status_code=500,
            detail=f"You cannot provide both 'name' and 'wid'",
        )
    if wid is not None or name is not None:
        # This needs to be explained in the API doc
        if wid:
            return dm.get_taxon_by_id(wid)
        return dm.get_taxa_with_name_matching(name)
        # TODO implement the rest

    return None


def references_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    structures = structures_from_structure_in_item(dm, item)

    if structures is None:
        return None

    return dm.get_references_of_structures(structures)


def references_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    taxa = taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return None

    return dm.get_references_of_taxa(taxa)


def structures_from_reference_in_item(dm: DataModel, item: Item) -> set[int] | None:
    references = references_from_reference_in_item(dm, item)

    if references is None:
        return None

    return dm.get_structures_of_references(references)


def structures_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    taxa = taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return None

    # TODO Set recursive=True to have all the structures from the parents too?
    #      We may have issues if we have a lot, and it will require a bit more work to get it with the db
    #      We could also have all the parenting relations in the DB and it would be much much faster
    out = set()
    for taxon in taxa:
        out.update(dm.get_structures_of_taxon(taxon))

    return out


def taxa_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    structures = structures_from_structure_in_item(dm, item)

    if structures is None:
        return None

    return dm.get_taxa_of_structures(structures)


def taxa_from_reference_in_item(dm: DataModel, item: Item) -> set[int] | None:
    references = references_from_reference_in_item(dm, item)

    if references is None:
        return None

    return dm.get_taxa_of_references(references)


def combine_and_filter_outputs(sets: list[set], limit: int) -> list[int]:
    non_none_outputs = [s for s in sets if s is not None]

    items = list(set.intersection(*non_none_outputs) if non_none_outputs else set())

    if limit == 0:
        return items
    else:
        return items[:limit]


def apply_limit(item: Item, items: list[Any] | set[Any]) -> list[Any]:
    lim = item.filter.limit

    if lim == 0:
        return list(items)
    else:
        return list(items)[:lim]


def get_triplets_for_item(item: Item, dm: DataModel) -> list[tuple[int, int, int]]:
    triplets_set = dm.get_triplets_for(
        reference_ids=references_from_reference_in_item(dm, item),
        structure_ids=structures_from_structure_in_item(dm, item),
        taxon_ids=taxa_from_taxon_in_item(dm, item),
    )

    return apply_limit(item, triplets_set)


def get_structures_for_item(item: Item, dm: DataModel) -> dict[int, str]:
    ids = combine_and_filter_outputs(
        [
            structures_from_structure_in_item(dm, item),
            structures_from_taxon_in_item(dm, item),
            structures_from_reference_in_item(dm, item),
        ],
        limit=item.filter.limit,
    )

    return dm.get_dict_of_sid_to_smiles(ids)


def get_taxa_for_item(item: Item, dm: DataModel) -> dict[int, str]:
    ids = combine_and_filter_outputs(
        [
            taxa_from_taxon_in_item(dm, item),
            taxa_from_structure_in_item(dm, item),
            taxa_from_reference_in_item(dm, item),
        ],
        limit=item.filter.limit,
    )

    return dm.get_dict_of_tid_to_taxon_name(ids)


def get_references_for_item(item: Item, dm: DataModel) -> dict[int, str]:
    ids = combine_and_filter_outputs(
        [
            references_from_reference_in_item(dm, item),
            references_from_structure_in_item(dm, item),
            references_from_taxon_in_item(dm, item),
        ],
        limit=item.filter.limit,
    )

    return dm.get_dict_of_rid_to_reference_doi(ids)
