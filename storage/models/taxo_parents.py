from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class TaxoParents(Base):
    __tablename__ = "taxo_parents"

    id: Mapped[int] = mapped_column(primary_key=True)
    parent_id: Mapped[int] = mapped_column(primary_key=True)
    distance: Mapped[int]

    __table_args__ = (
        Index("taxo_parent_id", "id"),
        Index("taxo_parent_parent_id", "parent_id"),
        Index("taxo_parent_ids", "id", "parent_id"),
        Index("taxo_parent_distance_with_id", "distance", "id"),
        Index("taxo_parent_distance_with_parent_id", "distance", "parent_id"),
    )

    def __repr__(self):
        return f"TaxoParents(id={self.id}, parent_id={self.parent_id}, distance={self.distance})"
