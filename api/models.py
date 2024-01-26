from typing import Dict, Optional, List, Union

from pydantic import BaseModel


# class ReferenceFilterItem(BaseModel):
#     TODO


class StructureFilterItem(BaseModel):
    substructure_search: bool = False
    similarity_level: float = 1.0


class TaxonFilterItem(BaseModel):
    taxon_children: bool = False


class FilterItem(BaseModel):
    # refrence: ReferenceFilterItem = ReferenceFilterItem()
    structure: StructureFilterItem = StructureFilterItem()
    taxon: TaxonFilterItem = TaxonFilterItem()
    limit: Optional[int] = None


class ReferenceItem(BaseModel):
    wid: Optional[int] = None
    doi: Optional[str] = None


class StructureItem(BaseModel):
    wid: Optional[int] = None
    molecule: Optional[str] = None


class TaxonItem(BaseModel):
    wid: Optional[int] = None
    name: Optional[str] = None


class Item(BaseModel):
    filter: FilterItem = FilterItem()
    reference: ReferenceItem = ReferenceItem()
    structure: StructureItem = StructureItem()
    taxon: TaxonItem = TaxonItem()
    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "reference": {
                        "wid": 44488598,
                        "doi": "10.1080/1057563021000040466",
                    },
                    "structure": {
                        "wid": 27151406,
                        "molecule": "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O",
                    },
                    "taxon": {"wid": 158572, "name": "Gentiana lutea"},
                    "filter": {
                        "structure": {
                            "substructure_search": True,
                            "similarity_level": 0.8,
                        },
                        "taxon": {"taxon_children": True},
                        "limit": 1000,
                    },
                }
            ]
        }
    }


class ReferenceObject(BaseModel):
    doi: str


class ReferenceResult(BaseModel):
    ids: List[int]
    objects: Optional[Dict[int, ReferenceObject]]
    count: Optional[int]
    description: Optional[str]


class StructureObject(BaseModel):
    smiles: str


class StructureResult(BaseModel):
    ids: List[int]
    objects: Optional[Dict[int, StructureObject]]
    count: Optional[int]
    description: Optional[str]


class TaxonObject(BaseModel):
    name: str


class TaxonResult(BaseModel):
    ids: List[int]
    objects: Optional[Dict[int, TaxonObject]]
    count: Optional[int]
    description: Optional[str]


class TripletResult(BaseModel):
    triplets: List[List[int]]
    references: Optional[Dict[int, ReferenceObject]]
    structures: Optional[Dict[int, StructureObject]]
    taxa: Optional[Dict[int, TaxonObject]]
    count: Optional[int]
    description: Optional[str]
