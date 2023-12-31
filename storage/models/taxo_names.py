from sqlalchemy import Index
from sqlalchemy.orm import Mapped, mapped_column

from storage.models.base import Base


class TaxoNames(Base):
    __tablename__ = "taxo_names"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str]

    __table_args__ = (
        Index("taxo_name_id", "id"),
        Index("taxo_name", "name"),
    )

    def __repr__(self):
        return f"TaxoNames(id={self.id}, name={self.name})"
