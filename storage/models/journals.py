from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class Journals(Base):
    __tablename__ = "journals"

    id: Mapped[int] = mapped_column(primary_key=True)
    title: Mapped[str]

    __table_args__ = (Index("journal_id", "id"),)

    def __repr__(self):
        return f"Journals(id={self.id}, title={self.title})"
