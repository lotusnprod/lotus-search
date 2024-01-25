from pydantic import BaseModel


class Item(BaseModel):
    structure_wid: int | None = None
    structure: str | None = None
    substructure_search: bool | None = None
    similarity_level: float = 1.0
    taxon_wid: int | None = None
    taxon_name: str | None = None
    reference_wid: int | None = None
    reference_doi: str | None = None
    limit: int = 1000
    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "structure_wid": 27151406,
                    "structure": "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O",
                    "substructure_search": True,
                    "similarity_level": 0.8,
                    "taxon_wid": 158572,
                    "taxon_name": "Gentiana lutea",
                    "reference_wid": 44488598,
                    "reference_doi": "10.1080/1057563021000040466",
                    "limit": 1000,
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
