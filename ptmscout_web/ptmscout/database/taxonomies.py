from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR
from ptmscout.database import Base, DBSession
from sqlalchemy.orm import relationship

class Taxonomy(Base):
    __tablename__='taxonomy'
    node_id = Column(Integer(10), primary_key=True)
    kingdom = Column(VARCHAR(1), index=True)
    name = Column(VARCHAR(100))
    strain = Column(VARCHAR(100))
    
    parent_id = Column(Integer(10), ForeignKey('taxonomy.node_id'))

    parent = relationship('Taxonomy', remote_side=[node_id])

    def __format_name(self):
        if self.strain:
            return "%s (%s)" % (self.name.strip(), self.strain.strip())
        return self.name.strip()

    formatted_name = property(__format_name)

    def save(self):
        DBSession.add(self)
        DBSession.flush()
    
class Species(Base):
    __tablename__='species'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    name = Column(VARCHAR(100), unique=True)
    taxon_id = Column(Integer(10), ForeignKey(Taxonomy.node_id))

    taxon = relationship(Taxonomy)

    def __init__(self, name):
        self.name = name

    def save(self):
        DBSession.add(self)
        DBSession.flush()
        
def getSpeciesByName(name):
    return DBSession.query(Species).filter_by(name=name).first()

def getAllSpecies():
    return DBSession.query(Species).all()

def getTaxonomyById(txid):
    return DBSession.query(Taxonomy).filter_by(node_id=txid).first()

def getTaxonByName(taxon, strain=None):
    return DBSession.query(Taxonomy).filter_by(name=taxon, strain=strain).first()
