from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class TaxoRankNames(Base):
    __tablename__ = "taxo_rank_names"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str]

    __table_args__ = (
        Index("taxo_rank_name_id", "id"),
        Index("taxo_rank_name", "name"),
    )

    def __repr__(self):
        return f"TaxoRankNames(id={self.id}, name={self.name})"
