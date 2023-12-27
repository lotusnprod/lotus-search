import logging

from fastapi import HTTPException

from api.models import Item
from model import DataModel

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def get_matching_references_from_reference_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns the WID of matching references."""
    references = None
    if item.reference_doi and item.reference_wid:
        raise HTTPException(
            status_code=500,
            detail=f"You cannot provide both 'reference_doi' and 'reference_wid'",
        )
    elif item.reference_wid is not None or item.reference_doi is not None:
        # This needs to be explained in the API doc
        if item.reference_wid:
            if item.reference_wid in dm.get_references():
                return {item.reference_wid}
            else:
                return set()
        else:
            if item.reference_doi:
                references = set(dm.get_references_with_doi(item.reference_doi))
                if references is None:
                    references = dm.get_references()
            else:
                references = dm.get_references()

        return references

    return references


def get_matching_structures_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns the WID of matching structures."""
    structures = None

    if item.structure and item.structure_wid:
        raise HTTPException(
            status_code=500,
            detail=f"You cannot provide both 'structure' and 'structure_wid'",
        )
    elif item.structure is not None or item.structure_wid is not None:
        # This needs to be explained in the API doc
        if item.structure_wid:
            if item.structure_wid in dm.structures_set():
                return {item.structure_wid}
            else:
                return set()
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

    return structures


def get_matching_taxa_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    """Returns the WID of matching taxa."""
    taxa = None
    if item.taxon_wid is not None or item.taxon_name is not None:
        # This needs to be explained in the API doc
        if item.taxon_wid:
            if item.taxon_wid in dm.get_taxa():
                return {item.taxon_wid}
            else:
                return set()
        else:
            if item.taxon_name:
                taxa = set(dm.get_taxa_with_name_containing(item.taxon_name))
                if taxa is None:
                    taxa = dm.get_taxa()
            else:
                taxa = dm.get_taxa()
        return taxa

    return taxa


# TODO WIP
# def get_matching_references_from_couple_in_item(dm: DataModel, item: Item) -> set[int]:
#     taxa = get_matching_taxa_from_taxon_in_item(dm, item)
#     # We need to get all the matching structures
#     structures = get_matching_structures_from_structure_in_item(dm, item)

#     if taxa is None:
#         return None

#     if structures is None:
#         return None

#     tax = set()
#     for taxon in taxa:
#         tax.update(dm.get_structures_of_taxon(taxon))

#     stru = set()
#     for structure in structures:
#         stru.update(dm.get_references_containing_structure(structure))

#     # TODO get triplets and intersect

#     return out


def get_matching_references_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    structures = get_matching_structures_from_structure_in_item(dm, item)

    if structures is None:
        return None

    return dm.get_references_of_structures(structures)


def get_matching_references_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return None

    return dm.get_references_of_taxa(taxa)


def get_matching_structures_from_reference_in_item(dm: DataModel, item: Item) -> set[int] | None:
    references = get_matching_references_from_reference_in_item(dm, item)

    if references is None:
        return None

    return dm.get_structures_of_references(references)


def get_matching_structures_from_taxon_in_item(dm: DataModel, item: Item) -> set[int] | None:
    taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return None

    # TODO Set recursive=True to have all the structures from the parents too?
    #      We may have issues if we have a lot, and it will require a bit more work to get it with the db
    #      We could also have all the parenting relations in the DB and it would be much much faster
    out = set()
    for taxon in taxa:
        out.update(dm.get_structures_of_taxon(taxon))

    return out


def get_matching_taxa_from_structure_in_item(dm: DataModel, item: Item) -> set[int] | None:
    structures = get_matching_structures_from_structure_in_item(dm, item)

    if structures is None:
        return None

    return dm.get_taxa_of_structures(structures)


def get_matching_taxa_from_reference_in_item(dm: DataModel, item: Item) -> set[int] | None:
    references = get_matching_references_from_reference_in_item(dm, item)

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
