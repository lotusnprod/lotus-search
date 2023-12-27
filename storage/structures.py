from sqlalchemy import Index, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column

from storage.base import Base


class Structures(Base):
    __tablename__ = "structures"

    id: Mapped[int] = mapped_column(primary_key=True)
    smiles: Mapped[str]

    __table_args__ = (
        Index("structure_id", "id"),
    )

    def __repr__(self):
        return f"Structures(id={self.id}, smiles={self.smiles})"
