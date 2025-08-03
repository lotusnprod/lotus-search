from sqlalchemy import Column, ForeignKey, Index, Integer

from storage.models.base import Base


class TaxoParents(Base):
    __tablename__ = "taxo_parents"

    id = Column(Integer, primary_key=True)
    child_id = Column(Integer, ForeignKey("triplets.taxon_id"))
    parent_id = Column(Integer)
    distance = Column(Integer)

    __table_args__ = (
        Index("taxo_parent_id", "id"),
        Index("taxo_parent_child_id", "child_id"),
        Index("taxo_parent_parent_id", "parent_id"),
        Index("taxo_parent_distance_id", "distance"),
    )

    def __repr__(self):
        return (
            f"TaxoParents(id={self.id}, child_id={self.child_id}, parent_id={self.parent_id}, distance={self.distance})"
        )
