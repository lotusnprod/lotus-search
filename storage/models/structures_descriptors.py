from sqlalchemy import Column, Float, ForeignKey, Index, Integer, String

from storage.models.base import Base


class StructuresDescriptors(Base):
    __tablename__ = "structures_descriptors"

    id = Column(Integer, primary_key=True)
    structure_id = Column(Integer, ForeignKey("structures.id"))
    descriptor_name = Column(String)
    descriptor_value = Column(Float)

    __table_args__ = (Index("descriptor_id", "descriptor_name"),)

    def __repr__(self):
        return f"StructuresDescriptors(id={self.id}, structure_id={self.structure_id}, descriptor_name={self.descriptor_name}, descriptor_value={self.descriptor_value})"
