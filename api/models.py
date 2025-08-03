from pydantic import BaseModel


class ReferenceOption(BaseModel):
    date_min: str | None = None
    date_max: str | None = None
    journal: str | None = None


class StructureOption(BaseModel):
    descriptors: dict[str, int] | None = None
    return_descriptors: bool = False
    substructure_search: bool = False
    similarity_level: float = 1.0
    sdf: bool = False


class TaxonOption(BaseModel):
    taxon_children: bool = False


class ReferenceItem(BaseModel):
    wid: int | None = None
    doi: str | None = None
    title: str | None = None
    option: ReferenceOption = ReferenceOption()
    # limit: Optional[int] = None


class StructureItem(BaseModel):
    wid: int | None = None
    molecule: str | None = None
    formula: str | None = None
    option: StructureOption = StructureOption()
    # limit: Optional[int] = None


class TaxonItem(BaseModel):
    wid: int | None = None
    name: str | None = None
    option: TaxonOption = TaxonOption()
    # limit: Optional[int] = None


class Item(BaseModel):
    reference: ReferenceItem = ReferenceItem()
    structure: StructureItem = StructureItem()
    taxon: TaxonItem = TaxonItem()
    limit: int | None = None
    modeEnum: str = "objects"
    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "reference": {
                        "wid": 44488598,
                        "doi": "10.1080/1057563021000040466",
                        "title": "Iridoids from Seeds of Gentiana Lutea",
                        "option": {
                            "date_min": "1999-12-31",
                            "date_max": "2024-01-01",
                            "journal": "Natural Product Research",
                        },
                    },
                    "structure": {
                        "wid": 27151406,
                        "molecule": "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O",
                        "formula": "C16H22O9",
                        "option": {
                            "descriptors": {"NumHAcceptors_min": 1},
                            "return_descriptors": False,
                            "substructure_search": True,
                            "similarity_level": 0.8,
                            "sdf": False,
                        },
                    },
                    "taxon": {
                        "wid": 158572,
                        "name": "Gentiana lutea",
                        "option": {"taxon_children": True},
                    },
                    "limit": 1000,
                    "modeEnum": "ids",
                },
            ],
        },
    }


class ReferenceObject(BaseModel):
    doi: str
    title: str
    date: str
    journal: str


class ReferenceResult(BaseModel):
    ids: list[int]
    objects: dict[int, ReferenceObject] | None = None
    count: int | None
    description: str | None


class StructureObject(BaseModel):
    smiles: str
    smiles_no_stereo: str
    inchi: str
    inchi_no_stereo: str
    inchikey: str
    inchikey_no_stereo: str
    formula: str
    descriptors: dict | None = None


class StructureResult(BaseModel):
    ids: list[int]
    objects: dict[int, StructureObject] | None = None
    sdf: str | None = None
    count: int | None
    description: str | None


class TaxonObject(BaseModel):
    name: str


class TaxonResult(BaseModel):
    ids: list[int]
    objects: dict[int, TaxonObject] | None = None
    count: int | None
    description: str | None


class TripletResult(BaseModel):
    triplets: list[list[int]]
    references: dict[int, ReferenceObject] | None = None
    structures: dict[int, StructureObject] | None = None
    taxa: dict[int, TaxonObject] | None = None
    count: int | None
    description: str | None


class AutocompleteTaxa(BaseModel):
    taxon_name: str


class DepictionStructure(BaseModel):
    structure: str
    highlight: str | None = None
