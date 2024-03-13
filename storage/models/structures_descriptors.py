from sqlalchemy import Column, Integer, String, Index, Float, ForeignKey
from sqlalchemy.orm import relationship

from storage.models.base import Base


class StructuresDescriptors(Base):
    __tablename__ = "structures_descriptors"

    # TODO this has not been tested, probably not working
    id = Column(Integer, primary_key=True)
    structure_id = Column(Integer, ForeignKey("structures.id"))
    structure = relationship("Structures", backref="descriptors")
    descriptor_name = Column(String)
    descriptor_value = Column(Float)

    __table_args__ = (Index("descriptor_id", "descriptor_name"),)

    def __repr__(self):
        return f"StructuresDescriptors(id={self.id}, structure_id={self.structure_id}, descriptor_name={self.descriptor_name}, descriptor_value={self.descriptor_value})"
