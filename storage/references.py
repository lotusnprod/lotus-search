from sqlalchemy import Index, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column

from storage.base import Base


class References(Base):
    __tablename__ = "references"

    id: Mapped[int] = mapped_column(primary_key=True)
    doi: Mapped[str]

    __table_args__ = (
        Index("reference_doi", "doi"),
        Index("reference_id", "id"),
    )

    def __repr__(self):
        return f"References(id={self.id}, doi={self.doi})"
