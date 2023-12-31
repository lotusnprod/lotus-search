from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class TaxoRanks(Base):
    __tablename__ = "taxo_ranks"

    id: Mapped[int] = mapped_column(primary_key=True)
    rank_id: Mapped[int] = mapped_column(primary_key=True)

    __table_args__ = (
        Index("taxo_rank_id", "id"),
        Index("taxo_rank_rank_id", "rank_id"),
    )

    def __repr__(self):
        return f"TaxoRanks(id={self.id}, rank_id={self.rank_id})"
