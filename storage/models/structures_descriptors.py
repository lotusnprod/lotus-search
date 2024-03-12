from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class StructuresDescriptors(Base):
    __tablename__ = "structures_descriptors"

    id: Mapped[int] = mapped_column(primary_key=True)
    smiles: Mapped[str]
    # TODO decide how to add the whole dict

    # TODO
    # __table_args__ = (Index("structure_id", "id"),)

    # TODO
    # def __repr__(self):
        # return f"StructuresDescriptors(id={self.id}, smiles={self.smiles})"
