from sqlalchemy import Column, ForeignKey, Index, Integer, String

from storage.models.base import Base


class TaxoRankNames(Base):
    __tablename__ = "taxo_rank_names"

    id = Column(Integer, ForeignKey("taxo_ranks.id"), primary_key=True)
    name = Column(String)

    __table_args__ = (
        Index("taxo_rank_name_id", "id"),
        Index("taxo_rank_name", "name"),
    )

    def __repr__(self):
        return f"TaxoRankNames(id={self.id}, name={self.name})"
