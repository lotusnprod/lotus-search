from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class Structures(Base):
    __tablename__ = "structures"

    id: Mapped[int] = mapped_column(primary_key=True)
    smiles: Mapped[str]
    smiles_no_stereo: Mapped[str]
    inchi: Mapped[str]
    inchi_no_stereo: Mapped[str]
    inchikey: Mapped[str]
    inchikey_no_stereo: Mapped[str]
    formula: Mapped[str]

    __table_args__ = (Index("structure_id", "id"),)

    def __repr__(self):
        return f"Structures(id={self.id}, smiles={self.smiles}, smiles_no_stereo={self.smiles_no_stereo}, inchi={self.inchi}, inchi_no_stereo={self.inchi_no_stereo}, inchikey={self.inchikey}, inchikey_no_stereo={self.inchikey_no_stereo}, formula={self.formula})"
