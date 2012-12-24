from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, UniqueConstraint, Table
from sqlalchemy.types import Integer, TEXT, VARCHAR, Enum, Text, Float, DateTime
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import null, or_, and_
from ptmscout.config import strings, settings
import taxonomies
import datetime

go_hierarchy_table = Table('GO_hierarchy', Base.metadata,
    Column('parent_id', Integer(10), ForeignKey('GO.id')),
    Column('child_id', Integer(10), ForeignKey('GO.id')))
    
expression_association_table = Table('protein_expression', Base.metadata,
    Column('id', Integer(10), primary_key=True, autoincrement=True),
    Column('protein_id', Integer(10), ForeignKey('protein.id')),
    Column('probeset_id', Integer(10), ForeignKey('expression.probeset_id')))


class GeneOntology(Base):
    __tablename__='GO'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    aspect = Column(Enum(['F','P','C']), default=null)
    GO = Column(VARCHAR(10))
    term = Column(Text)
    date = Column(DateTime)
    version = Column(VARCHAR(10))
    UniqueConstraint('aspect', 'GO', name="uniqueEntry")
    
    
    children = relationship("GeneOntology", secondary=go_hierarchy_table, backref="parents",
                        primaryjoin=id==go_hierarchy_table.c.parent_id,
                        secondaryjoin=id==go_hierarchy_table.c.child_id)
    
    def __init__(self):
        self.date = datetime.datetime.now()
    
#    def getURL(self):
#        return settings.accession_urls['GO'] % (self.GO)

    def save(self):
        DBSession.add(self)

class GeneOntologyEntry(Base):
    __tablename__='protein_GO'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    GO_id = Column(Integer(10), ForeignKey('GO.id'))
    date = Column(DateTime)
    
    GO_term = relationship("GeneOntology")

    def save(self):
        DBSession.add(self)

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
    
    def save(self):
        DBSession.add(self)

class Protein(Base):
    __tablename__='protein'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    sequence = Column(TEXT)
    acc_gene = Column(VARCHAR(30))
    name = Column(VARCHAR(100))
    date = Column(DateTime)
    species_id = Column(Integer(10), ForeignKey('species.id'))
    
    accessions = relationship("ProteinAccession", order_by=ProteinAccession.type)
    domains = relationship("ProteinDomain")
    
    species = relationship("Species")
    GO_terms = relationship("GeneOntologyEntry")
    expression_probes = relationship("ExpressionProbeset", secondary=expression_association_table)
    
    def __init__(self):
        self.date = datetime.datetime.now()
    
    def saveProtein(self):
        DBSession.add(self)
        DBSession.flush()
        
    def hasAccession(self, acc):
        lower_acc = acc.lower()
        for dbacc in self.accessions:
            if dbacc.value.lower() == lower_acc:
                return True
        return False
        
    def addTaxonomy(self, taxon):
        if taxon not in self.taxonomy:
            self.taxonomy.append(taxon)
    
    def addGoTerm(self, GO_term, date_added=datetime.datetime.now()):
        goe = GeneOntologyEntry()
        goe.GO_term = GO_term
        goe.version = date_added
        self.GO_terms.append(goe)
    
    def getGOIds(self):
        return set([ goe.GO_term.GO for goe in self.GO_terms ])
     
    def hasGoTerm(self, GO_id):
        for goe in self.GO_terms:
            if goe.GO_term.GO == GO_id:
                return True
        
        return False

class NoSuchProtein(Exception):
    def __init__(self, prot):
        self.prot = prot
    
    def __str__(self):
        return "No such protein: %s" % (str(self.prot))

def getGoAnnotationById(goId):
    value = DBSession.query(GeneOntology).filter_by(GO=goId).first()
    
    return value

def getProteinById(pid):
    value = DBSession.query(Protein).filter_by(id=pid).first()
    
    if value == None:
        raise NoSuchProtein(pid)
    
    return value

def getProteinsByAccession(accessions, species=None):
    if species == None:
        q = DBSession.query(Protein).join(Protein.accessions).filter(or_(Protein.acc_gene.in_(accessions), ProteinAccession.value.in_(accessions)))
    else:
        q = DBSession.query(Protein).join(Protein.accessions).join(Protein.species).filter(or_(Protein.acc_gene.in_(accessions), ProteinAccession.value.in_(accessions)), taxonomies.Species.name == species)
    
    return q.all()

def searchProteins(search, species=None, page=None):
    search = "%" + search + "%"

    q = DBSession.query(Protein.id).join(Protein.accessions).join(Protein.species)
    or_clause = or_(Protein.acc_gene.like(search),
            ProteinAccession.value.like(search),
            Protein.name.like(search),
            Protein.sequence.like(search))

    if species:
        q = q.filter(and_(or_clause, taxonomies.Species.name==species)).distinct()
    else:
        q = q.filter(or_clause).distinct()

    if page==None:
        sq = q.subquery()
    else:
        limit, offset = page
        sq = q.limit(limit).offset(offset).subquery()

    fq = DBSession.query(Protein).join(sq, Protein.id == sq.c.id)
    return fq.count(), fq.all()

def getProteinDomain(prot_id, pep_site):
    return DBSession.query(ProteinDomain).filter(and_(ProteinDomain.protein_id==prot_id, ProteinDomain.start <= pep_site, pep_site <= ProteinDomain.stop)).first()    


def getProteinBySequence(seq, species):
    return DBSession.query(Protein).join(Protein.species).filter(Protein.sequence==seq, taxonomies.Species.name==species).first()


