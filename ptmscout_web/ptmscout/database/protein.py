from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, UniqueConstraint, Table
from sqlalchemy.types import Integer, TEXT, VARCHAR, Enum, Text, Float
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import null, or_
from ptmscout import strings, config

go_association_table = Table('protein_GO', Base.metadata,
    Column('protein_id', Integer(10), ForeignKey('protein.id')),
    Column('GO_id', Integer(10), ForeignKey('GO.id')),
    Column('version', VARCHAR(10))
)

class Species(Base):
    __tablename__='species'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    name = Column(VARCHAR(100), unique=True)
    
    def __init__(self, name):
        self.name = name

class GeneOntology(Base):
    __tablename__='GO'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    aspect = Column('aspect', Enum(['F','P','C']), default=null)
    GO = Column(VARCHAR(10))
    term = Column(Text)
    version = Column(VARCHAR(10), default='0')
    UniqueConstraint('aspect', 'GO', name="uniqueEntry")
    
#    def getURL(self):
#        return config.accession_urls['GO'] % (self.GO)

class Accession(Base):
    __tablename__='acc'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    type = Column(VARCHAR(30))
    value = Column(VARCHAR(45))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
    def getType(self):
        return strings.accession_type_strings[self.type]
    
    def getURL(self):
        if self.type in config.accession_urls:
            return config.accession_urls[self.type] % (self.value)
        return None
        
    def getAccessionName(self):
        return self.value
    
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
    acc_gene = Column(VARCHAR(30))
    name = Column(VARCHAR(100))
    date = Column(VARCHAR(7))
    species_id = Column(Integer(10), ForeignKey('species.id'))
    
    accessions = relationship("Accession")
    GO_terms = relationship("GeneOntology", secondary=go_association_table)
    domains = relationship("Domain")
    species = relationship("Species")
    
    def __init__(self):
        pass
    
    def saveProtein(self):
        DBSession.add(self)

class NoSuchProtein(Exception):
    def __init__(self, prot):
        self.prot = prot
    
    def __str__(self):
        return "No such protein: %s" % (str(self.prot))

def getProteinById(pid):
    value = DBSession.query(Protein).filter_by(id=pid).first()
    
    if value == None:
        raise NoSuchProtein(pid)
    
    return value

def getProteinsByAccession(accessions, species=None):
    if species == None:
        q = DBSession.query(Protein).join(Protein.accessions).filter(or_(Protein.acc_gene.in_(accessions), Accession.value.in_(accessions)))
    else:
        q = DBSession.query(Protein).join(Protein.accessions).join(Protein.species).filter(or_(Protein.acc_gene.in_(accessions), Accession.value.in_(accessions)), Species.name == species)
    
    return q.all()

def getSpeciesByName(name):
    return DBSession.query(Species).filter_by(name=name).first()

def getAllSpecies():
    return DBSession.query(Species).all()