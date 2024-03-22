from sqlalchemy import Column, ForeignKey, Index, Integer, String

from storage.models.base import Base


class Journals(Base):
    __tablename__ = "journals"

    id = Column(Integer, ForeignKey("references.journal"), primary_key=True)
    title = Column(String)

    __table_args__ = (Index("journal_id", "id"),)

    def __repr__(self):
        return f"Journals(id={self.id}, title={self.title})"
