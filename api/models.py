from typing import Optional, Union

from pydantic import BaseModel


class Reference(BaseModel):
    wid: Optional[int] = None
    doi: Optional[str] = None


class Structure(BaseModel):
    wid: Optional[int] = None
    molecule: Optional[str] = None


class Taxon(BaseModel):
    wid: Optional[int] = None
    name: Optional[str] = None


class StructuralFilter(BaseModel):
    substructure_search: bool = False
    similarity_level: float = 1.0


class TaxalFilter(BaseModel):
    taxon_children: bool = False


class Filter(BaseModel):
    structural: StructuralFilter = StructuralFilter()
    taxal: TaxalFilter = TaxalFilter()
    limit: Optional[int] = None


class Item(BaseModel):
    reference: Reference = Reference()
    structure: Structure = Structure()
    taxon: Taxon = Taxon()
    filter: Filter = Filter()

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
                        "structural": {
                            "substructure_search": True,
                            "similarity_level": 0.8,
                        },
                        "taxal": {"taxon_children": True},
                        "limit": 1000,
                    },
                }
            ]
        }
    }


class ReferenceInfo(BaseModel):
    doi: str


class ReferenceResult(BaseModel):
    ids: list[int]
    references: dict[int, ReferenceInfo]
    count: int
    description: str


class StructureInfo(BaseModel):
    smiles: str


class StructureResult(BaseModel):
    ids: list[int]
    structures: dict[int, StructureInfo]
    count: int
    description: str


class TaxonInfo(BaseModel):
    name: str


class TaxonResult(BaseModel):
    ids: list[int]
    taxa: dict[int, TaxonInfo]
    count: int
    description: str


class TripletResult(BaseModel):
    triplets: list[list[int]]
    references: dict[int, ReferenceInfo]
    structures: dict[int, StructureInfo]
    taxa: dict[int, TaxonInfo]
    count: int
    description: str


class AutocompleteTaxa(BaseModel):
    taxon_name: str


class StructureDepiction(BaseModel):
    structure: str
    highlight: str | None = None
