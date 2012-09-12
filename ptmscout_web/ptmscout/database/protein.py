from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, UniqueConstraint, Table
from sqlalchemy.types import Integer, TEXT, VARCHAR, Enum, Text, Float
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import null

go_association_table = Table('protein_GO', Base.metadata,
    Column('protein_id', Integer(10), ForeignKey('protein.id')),
    Column('GO_id', Integer(10), ForeignKey('GO.id')),
    Column('version', VARCHAR(10))
)

class GeneOntology(Base):
    __tablename__='GO'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    aspect = Column('aspect', Enum(['F','P','C']), default=null)
    GO = Column(VARCHAR(10))
    term = Column(Text)
    version = Column(VARCHAR(10), default='0')
    UniqueConstraint('aspect', 'GO', name="uniqueEntry")

class Accession(Base):
    __tablename__='acc'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    type = Column(VARCHAR(30))
    value = Column(VARCHAR(45))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
class Domain(Base):    
    __tablename__='domain'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    label = Column(VARCHAR(45))
    start = Column(Integer(10))
    stop = Column(Integer(10))
    p_value = Column(Float)
    source = Column(Text, default='pfam')
    params = Column(VARCHAR(45))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    version = Column(Integer(11))

class Protein(Base):
    __tablename__='protein'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    sequence = Column(TEXT)
    species = Column(VARCHAR(100))
    acc_gene = Column(VARCHAR(30))
    name = Column(VARCHAR(100))
    date = Column(VARCHAR(7))
    
    accessions = relationship("Accession")
    GO_terms = relationship("GeneOntology", secondary=go_association_table)
    domains = relationship("Domain")
    
    def __init__(self):
        print "Got protein: ", id
        pass
    
    def saveProtein(self):
        DBSession.add(self)

class NoSuchProtein(Exception):
    def __init__(self, pid):
        self.pid = pid
    
    def __str__(self):
        return "No such protein: %d" % (self.pid)

def getProteinById(pid):
    value = DBSession.query(Protein).filter_by(id=pid).first()
    
    if value == None:
        raise NoSuchProtein(pid)
    
    return value