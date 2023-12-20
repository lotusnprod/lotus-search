import logging

from api.models import (CoupleResult, Item, ReferenceInfo, ReferenceResult,
                        StructureInfo, StructureResult, TaxonInfo, TaxonResult)
from model import DataModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def get_matching_references_from_reference_in_item(
    dm: DataModel, item: Item
) -> set[int] | None:
    """Returns the WID of matching references."""
    if item.reference_wid is None and item.reference_doi is None:
        return None
    else:
        # This needs to be explained in the API doc
        if item.reference_wid:
            if item.reference_wid in dm.get_references():
                return {item.reference_wid}
        else:
            if item.reference_doi:
                references = set(dm.get_references_with_doi(item.reference_doi))
                if references is None:
                    references = dm.get_references()
            else:
                references = dm.get_references()

        return references


def get_matching_structures_from_structure_in_item(
    dm: DataModel, item: Item
) -> set[int]:
    """Returns the WID of matching structures."""
    if item.structure is None and item.structure_wid is None:
        return None
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
    """Returns the WID of matching taxa."""
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


# TODO WIP
# def get_matching_references_from_couple_in_item(dm: DataModel, item: Item) -> set[int]:
#     # We need to get all the matching taxa
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

#     # TODO get couples and intersect

#     return out


def get_matching_references_from_structure_in_item(
    dm: DataModel, item: Item
) -> set[int]:
    # We need to get all the matching structures
    structures = get_matching_structures_from_structure_in_item(dm, item)

    if structures is None:
        return set()

    out = set()
    for structure in structures:
        out.update(dm.get_references_containing_structure(structure))

    return out


def get_matching_references_from_taxon_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching taxa
    taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return set()

    out = set()
    for taxon in taxa:
        out.update(dm.get_references_containing_taxon(taxon))

    return out


def get_matching_structures_from_reference_in_item(
    dm: DataModel, item: Item
) -> set[int]:
    # We need to get all the matching references
    references = get_matching_references_from_reference_in_item(dm, item)

    if references is None:
        return set()

    out = set()
    for reference in references:
        out.update(dm.get_structures_of_reference(reference))

    return out


def get_matching_structures_from_taxon_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching taxa
    taxa = get_matching_taxa_from_taxon_in_item(dm, item)

    if taxa is None:
        return set()

    # Set recursive=True to have all the structures from the parents too
    out = set()
    for taxon in taxa:
        out.update(dm.get_structures_of_taxon(taxon))

    return out


def get_matching_taxa_from_structure_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching structures
    structures = get_matching_structures_from_structure_in_item(dm, item)

    if structures is None:
        return set()

    out = set()
    for structure in structures:
        out.update(dm.get_taxa_containing_structure(structure))

    return out


def get_matching_taxa_from_reference_in_item(dm: DataModel, item: Item) -> set[int]:
    # We need to get all the matching references
    references = get_matching_references_from_reference_in_item(dm, item)

    if references is None:
        return set()

    out = set()
    for reference in references:
        out.update(dm.get_taxa_of_reference(references))

    return out
