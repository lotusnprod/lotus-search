from typing import Dict, List, Optional, Union

from pydantic import BaseModel

# class ReferenceOption(BaseModel):
#     TODO


class StructureOption(BaseModel):
    substructure_search: bool = False
    similarity_level: float = 1.0


class TaxonOption(BaseModel):
    taxon_children: bool = False


class ReferenceItem(BaseModel):
    wid: Optional[int] = None
    doi: Optional[str] = None
    # option: ReferenceOption = ReferenceOption()
    # limit: Optional[int] = None


class StructureItem(BaseModel):
    wid: Optional[int] = None
    molecule: Optional[str] = None
    option: StructureOption = StructureOption()
    # limit: Optional[int] = None


class TaxonItem(BaseModel):
    wid: Optional[int] = None
    name: Optional[str] = None
    option: TaxonOption = TaxonOption()
    # limit: Optional[int] = None


class Item(BaseModel):
    reference: ReferenceItem = ReferenceItem()
    structure: StructureItem = StructureItem()
    taxon: TaxonItem = TaxonItem()
    limit: Optional[int] = None
    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "reference": {
                        "wid": 44488598,
                        "doi": "10.1080/1057563021000040466",
                        # "option": {
                        #     "TODO": True,
                        # },
                        # "limit": 1000
                    },
                    "structure": {
                        "wid": 27151406,
                        "molecule": "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O",
                        "option": {
                            "substructure_search": True,
                            "similarity_level": 0.8,
                        },
                        # "limit": 1000
                    },
                    "taxon": {
                        "wid": 158572,
                        "name": "Gentiana lutea",
                        "option": {"taxon_children": True},
                        # "limit": 1000
                    },
                    "limit": 1000,
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


class AutocompleteTaxa(BaseModel):
    taxon_name: str


class StructureDepiction(BaseModel):
    structure: str
    highlight: str | None = None
