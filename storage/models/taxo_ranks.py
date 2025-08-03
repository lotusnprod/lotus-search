from sqlalchemy import Column, ForeignKey, Index, Integer

from storage.models.base import Base


class TaxoRanks(Base):
    __tablename__ = "taxo_ranks"

    id = Column(Integer, ForeignKey("triplets.taxon_id"), primary_key=True)
    rank_id = Column(Integer)

    __table_args__ = (
        Index("taxo_rank_id", "id"),
        Index("taxo_rank_rank_id", "rank_id"),
    )

    def __repr__(self):
        return f"TaxoRanks(id={self.id}, rank_id={self.rank_id})"
