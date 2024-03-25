from sqlalchemy import Column, ForeignKey, Index, Integer, String

from storage.models.base import Base


class References(Base):
    __tablename__ = "references"

    id = Column(Integer, ForeignKey("triplets.reference_id"), primary_key=True)
    doi = Column(String)
    title = Column(String)
    date = Column(String)
    journal = Column(Integer)

    __table_args__ = (
        Index("reference_doi", "doi"),
        Index("reference_id", "id"),
        Index("reference_title", "title"),
        Index("reference_date", "date"),
        Index("reference_journal", "journal"),
    )

    def __repr__(self):
        return f"References(id={self.id}, doi={self.doi}, title={self.title}, date={self.date}, journal={self.journal})"
