from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class References(Base):
    __tablename__ = "references"

    id: Mapped[int] = mapped_column(primary_key=True)
    doi: Mapped[str]
    # TODO uncomment when ready
    # title: Mapped[str]
    # date: Mapped[str]
    # journal: Mapped[int]

    __table_args__ = (
        Index("reference_doi", "doi"),
        Index("reference_id", "id"),
        # Index("reference_title", "title"),
        # Index("reference_date", "date"),
        # Index("reference_journal", "journal"),
    )

    def __repr__(self):
        return f"References(id={self.id}, doi={self.doi})"
        # return f"References(id={self.id}, doi={self.doi}, title={self.title}, date={self.date}, journal={self.journal})"
