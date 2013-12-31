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
    
    def hasChild(self, goId):
        goId = goId.lower()
        for c in self.children:
            if c.GO.lower() == goId:
                return True
        return False

    def fullName(self):
        return "%s: %s" % (self.GO, self.term)
#    def getURL(self):
#        return settings.accession_urls['GO'] % (self.GO)

    def save(self):
        DBSession.add(self)
        DBSession.flush()

class GeneOntologyEntry(Base):
    __tablename__='protein_GO'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    GO_id = Column(Integer(10), ForeignKey('GO.id'))
    date = Column(DateTime)
    
    GO_term = relationship("GeneOntology", lazy='joined')

class ProteinScansite(Base):
    __tablename__='protein_scansite'
    id = Column(Integer(10), autoincrement=True, primary_key=True)
    
    source = Column(VARCHAR(40), default='scansite')
    value = Column(VARCHAR(20))
    score = Column(Float)
    percentile = Column(Float)
    site_pos = Column(Integer(10))
    
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
    UniqueConstraint('source', 'value', 'peptide_id', name="UNIQUE_pepId_source_value")

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

    def hasSite(self, site_pos):
        return self.start <= site_pos and site_pos <= self.stop


class ProteinRegion(Base):
    __tablename__='protein_regions'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    type = Column(VARCHAR(20))
    label = Column(VARCHAR(100))
    source = Column(Enum(['predicted', 'parsed', 'uniprot', 'ncbi']))
    start = Column(Integer(10))
    stop = Column(Integer(10))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))

    def __init__(self, type, label, source, start, stop, protein_id=None):
        self.type = type
        self.label = label
        self.source = source
        self.start = start
        self.stop = stop
        self.protein_id = protein_id

    def __eq__(self, o):
        c0 = self.type == o.type
        c1 = self.label == o.label
        c2 = self.start == o.start
        c3 = self.stop == o.stop
        return c0 and c1 and c2 and c3

    def hasSite(self, site_pos):
        return self.start <= site_pos and site_pos <= self.stop

class Protein(Base):
    __tablename__='protein'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    sequence = Column(TEXT)
    acc_gene = Column(VARCHAR(30))
    locus = Column(VARCHAR(30))
    name = Column(VARCHAR(100))
    date = Column(DateTime)
    species_id = Column(Integer(10), ForeignKey('species.id'))
    
    accessions = relationship("ProteinAccession", order_by=ProteinAccession.type, cascade="all,delete-orphan")
    domains = relationship("ProteinDomain")
    
    species = relationship("Species")
    GO_terms = relationship("GeneOntologyEntry", lazy='joined')
    expression_probes = relationship("ExpressionProbeset", secondary=expression_association_table)
    mutations = relationship("Mutation", cascade="all,delete-orphan")
    regions = relationship("ProteinRegion")
    scansite = relationship("ProteinScansite")

    def __init__(self):
        self.date = datetime.datetime.now()
    
    def hasPrediction(self, source, value, site_pos):
        for pred in self.scansite:
            if source == pred.source and \
                value == pred.value and \
                site_pos == pred.site_pos:
                return True
        return False
    
    def saveProtein(self):
        DBSession.add(self)
        DBSession.flush()
   
    def saveNoFlush(self):
        DBSession.add(self)

    def getGeneName(self):
        if self.acc_gene != None and self.acc_gene != '':
            return self.acc_gene
        return self.locus

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
        goe.date = date_added

        self.GO_terms.append(goe)
    
    def getGOIds(self):
        return set([ goe.GO_term.GO for goe in self.GO_terms ])
     
    def hasGoTerm(self, GO_id):
        for goe in self.GO_terms:
            if goe.GO_term.GO.lower() == GO_id.lower():
                return True
        return False

    def hasMutation(self, m):
        compare = [m.equals(m2) for m2 in self.mutations]
        return reduce(bool.__or__, compare, False)

    def hasRegion(self, region):
        compare = [region == r2 for r2 in self.regions]
        return reduce(bool.__or__, compare, False)

    def get_kmer(self, site_pos, k=7):
        prot_seq = self.sequence
        site_pos = site_pos - 1
        low_bound = max([site_pos-k, 0])
        high_bound = min([len(prot_seq), site_pos+k+1])

        pep_aligned = prot_seq[low_bound:site_pos] + prot_seq[site_pos].lower() + prot_seq[site_pos+1:high_bound]
        
        if site_pos-k < 0:
            pep_aligned = (" " * (k - site_pos)) + pep_aligned
        if site_pos+k+1 > len(prot_seq):
            pep_aligned = pep_aligned + (" " * (site_pos + k+1 - len(prot_seq)))

        return pep_aligned

    def get_domain(self, site_pos):
        for d in self.domains:
            if d.start <= site_pos and site_pos < d.stop:
                return d

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

def getProteinsByGene(gene_name, species=None):
    if species == None:
        q = DBSession.query(Protein).join(Protein.accessions).filter( or_( Protein.acc_gene==gene_name, ProteinAccession.value==gene_name ) )
    else:
        q = DBSession.query(Protein).join(Protein.accessions).join(Protein.species).filter( or_( Protein.acc_gene==gene_name, ProteinAccession.value==gene_name ), taxonomies.Species.name == species)
    
    return q.all()

def searchProteins(search=None, species=None, sequence=None, page=None, exp_id=None, includeNames=False):
    q = DBSession.query(Protein.id).join(Protein.accessions).join(Protein.species)

    clause = "1=1"
    if search:
        search = ( "%" + search + "%" ) if search else "%"
        if includeNames:
            clause = or_(Protein.acc_gene.like(search),
                    ProteinAccession.value.like(search),
                    Protein.name.like(search))
        else:
            clause = or_(Protein.acc_gene.like(search),
                    ProteinAccession.value.like(search))

    if sequence:
        clause = and_(clause, Protein.sequence.op('regexp')(sequence))

    if exp_id:
        from ptmscout.database import modifications
        sq = modifications.queryProteinsByExperiment(exp_id).subquery()
        q = q.join(sq, Protein.id == sq.c.protein_id)

    if species:
        clause = and_(clause, taxonomies.Species.name==species)

    q = q.filter(clause).distinct().order_by(Protein.acc_gene)

    if page==None:
        sq = q.subquery()
    else:
        limit, offset = page
        sq = q.limit(limit).offset(offset).subquery()

    fq = DBSession.query(Protein).join(sq, Protein.id == sq.c.id)
    return q.count(), fq.all()

def getProteinsByExperiment(exp_id, page=None):
    from ptmscout.database import modifications

    sq = modifications.queryProteinsByExperiment(exp_id).subquery()
    q = DBSession.query(Protein.id).join(sq, Protein.id == sq.c.protein_id).order_by(Protein.acc_gene)
    result_size = q.count()

    if page:
        limit, offset = page
        q = q.limit(limit).offset(offset)
    sq = q.subquery()

    return result_size, DBSession.query(Protein).join(sq, Protein.id == sq.c.id).all()

def getProteinDomain(prot_id, site_pos):
    return DBSession.query(ProteinDomain).filter(and_(ProteinDomain.protein_id==prot_id, ProteinDomain.start <= site_pos, site_pos <= ProteinDomain.stop)).first()    


def getProteinBySequence(seq, species):
    return DBSession.query(Protein).join(Protein.species).filter(Protein.sequence==seq, taxonomies.Species.name==species).first()

def getProteinsByAccession(accession):
    return DBSession.query(Protein).join(ProteinAccession).filter(ProteinAccession.value == accession).all()

def getAllProteins():
    return DBSession.query(Protein).all()

def countProteins():
    return DBSession.query(Protein).count()
