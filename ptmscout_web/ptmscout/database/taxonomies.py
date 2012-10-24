from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR
from ptmscout.database import Base

class Taxonomy(Base):
    __tablename__='taxonomy'
    node_id = Column(Integer(10), primary_key=True)
    kingdom = Column(VARCHAR(1), index=True)
    name = Column(VARCHAR(100))
    strain = Column(VARCHAR(100))
    
class Species(Base):
    __tablename__='species'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    name = Column(VARCHAR(100), unique=True)
    taxon_id = Column(Integer(10), ForeignKey(Taxonomy.node_id))
    
    def __init__(self, name):
        self.name = name