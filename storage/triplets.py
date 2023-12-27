from sqlalchemy import Index, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column

from storage.base import Base


class Triplets(Base):
    __tablename__ = "triplets"

    id: Mapped[int] = mapped_column(primary_key=True)
    reference_id: Mapped[int]
    structure_id: Mapped[int]
    taxon_id: Mapped[int]

    __table_args__ = (
        Index("r", "reference_id"),
        Index("s", "structure_id"),
        Index("t", "taxon_id"),
        Index("rs", "reference_id", "structure_id"),
        Index("rt", "reference_id", "taxon_id"),
        Index("st", "structure_id", "taxon_id"),
        Index("rst", "reference_id", "structure_id", "taxon_id"),
        UniqueConstraint('reference_id', 'structure_id', 'taxon_id', name='unique_rst')
    )

    def __repr__(self):
        return f"Triplets(id={self.id}, reference={self.reference_id}, structure={self.structure_id}, taxon={self.taxon_id})"
