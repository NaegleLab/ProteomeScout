from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, UniqueConstraint, Table
from sqlalchemy.types import Integer, TEXT, VARCHAR, Enum, Text, Float, DATETIME
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import null, or_
from ptmscout.config import strings, settings
from ptmscout.database.taxonomies import Species
import datetime

go_association_table = Table('protein_GO', Base.metadata,
    Column('protein_id', Integer(10), ForeignKey('protein.id')),
    Column('GO_id', Integer(10), ForeignKey('GO.id')),
    Column('version', VARCHAR(10)))

expression_association_table = Table('protein_expression', Base.metadata,
    Column('id', Integer(10), primary_key=True, autoincrement=True),
    Column('protein_id', Integer(10), ForeignKey('protein.id')),
    Column('probeset_id', Integer(10), ForeignKey('expression_ann.probeset_id')))

class GeneOntology(Base):
    __tablename__='GO'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    aspect = Column(Enum(['F','P','C']), default=null)
    GO = Column(VARCHAR(10))
    term = Column(Text)
    version = Column(VARCHAR(10), default='0')
    UniqueConstraint('aspect', 'GO', name="uniqueEntry")
    
#    def getURL(self):
#        return settings.accession_urls['GO'] % (self.GO)

class ProteinAccession(Base):
    __tablename__='protein_acc'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    type = Column(VARCHAR(30))
    value = Column(VARCHAR(45))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
    def getType(self):
        return strings.accession_type_strings[self.type]
    
    def getURL(self):
        if self.type in settings.accession_urls:
            return settings.accession_urls[self.type] % (self.value)
        return None
        
    def getAccessionName(self):
        return self.value
    
class ProteinDomain(Base):    
    __tablename__='protein_domain'
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
    date = Column(DATETIME)
    species_id = Column(Integer(10), ForeignKey('species.id'))
    
    accessions = relationship("ProteinAccession", order_by=ProteinAccession.type)
    domains = relationship("ProteinDomain")
    
    species = relationship("Species")
    GO_terms = relationship("GeneOntology", secondary=go_association_table)
    expression_probes = relationship("ExpressionProbeset", secondary=expression_association_table)
    
    def __init__(self):
        self.date = datetime.datetime.now()
    
    def saveProtein(self):
        DBSession.add(self)
        DBSession.flush()
        
    def addTaxonomy(self, taxon):
        if taxon not in self.taxonomy:
            self.taxonomy.append(taxon)
        

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
        q = DBSession.query(Protein).join(Protein.accessions).filter(or_(Protein.acc_gene.in_(accessions), ProteinAccession.value.in_(accessions)))
    else:
        q = DBSession.query(Protein).join(Protein.accessions).join(Protein.species).filter(or_(Protein.acc_gene.in_(accessions), ProteinAccession.value.in_(accessions)), Species.name == species)
    
    return q.all()
