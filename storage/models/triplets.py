from sqlalchemy import Column, Index, Integer, UniqueConstraint

from storage.models.base import Base


class Triplets(Base):
    __tablename__ = "triplets"

    id = Column(Integer, primary_key=True)
    reference_id = Column(Integer)
    structure_id = Column(Integer)
    taxon_id = Column(Integer)

    __table_args__ = (
        Index("r", "reference_id"),
        Index("s", "structure_id"),
        Index("t", "taxon_id"),
        Index("rs", "reference_id", "structure_id"),
        Index("rt", "reference_id", "taxon_id"),
        Index("st", "structure_id", "taxon_id"),
        Index("rst", "reference_id", "structure_id", "taxon_id"),
        UniqueConstraint("reference_id", "structure_id", "taxon_id", name="unique_rst"),
    )

    def __repr__(self):
        return f"Triplets(id={self.id}, reference={self.reference_id}, structure={self.structure_id}, taxon={self.taxon_id})"
