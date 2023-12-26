from sqlalchemy.orm import Mapped, mapped_column

from storage.base import Base


class SchemaVersion(Base):
    __tablename__ = "schema_version"

    version: Mapped[int] = mapped_column(primary_key=True)

    def __repr__(self):
        return f"SchemaVersion(version={self.version})"
