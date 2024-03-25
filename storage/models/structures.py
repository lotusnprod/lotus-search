from sqlalchemy import Column, ForeignKey, Index, Integer, String

from storage.models.base import Base


class Structures(Base):
    __tablename__ = "structures"

    id = Column(Integer, ForeignKey("triplets.structure_id"), primary_key=True)
    smiles = Column(String)
    smiles_no_stereo = Column(String)
    inchi = Column(String)
    inchi_no_stereo = Column(String)
    inchikey = Column(String)
    inchikey_no_stereo = Column(String)
    formula = Column(String)

    __table_args__ = (Index("structure_id", "id"),)

    def __repr__(self):
        return f"Structures(id={self.id}, smiles={self.smiles}, smiles_no_stereo={self.smiles_no_stereo}, inchi={self.inchi}, inchi_no_stereo={self.inchi_no_stereo}, inchikey={self.inchikey}, inchikey_no_stereo={self.inchikey_no_stereo}, formula={self.formula})"
