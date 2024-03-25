from sqlalchemy import Column, ForeignKey, Index, Integer, String

from storage.models.base import Base


class TaxoNames(Base):
    __tablename__ = "taxo_names"

    id = Column(Integer, ForeignKey("triplets.taxon_id"), primary_key=True)
    name = Column(String)

    __table_args__ = (
        Index("taxo_name_id", "id"),
        Index("taxo_name", "name"),
    )

    def __repr__(self):
        return f"TaxoNames(id={self.id}, name={self.name})"
