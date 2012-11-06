from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, Table, UniqueConstraint
from sqlalchemy.types import Integer, VARCHAR, CHAR, Float, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import and_, or_

MS_phosphopep = Table('MS_phosphopep', Base.metadata,
                      Column('id', Integer(10), primary_key=True),
                      Column('MS_id', Integer(10), ForeignKey('MS.id')),
                      Column('phosphopep_id', Integer(10), ForeignKey('phosphopep.id')))

PTM_taxon = Table('PTM_taxonomy', Base.metadata,
                    Column('PTM_id', Integer(10), ForeignKey('PTM.id')),
                    Column('taxon_id', Integer(10), ForeignKey('taxonomy.node_id')))

class PTMkeyword(Base):
    __tablename__ = 'PTM_keywords'
    id = Column(Integer(10), autoincrement=True, primary_key=True)
    PTM_id = Column(Integer(10), ForeignKey('PTM.id'))
    keyword = Column(VARCHAR(100))
    
class PTM(Base):
    __tablename__ = 'PTM'
    
    id = Column(Integer(10), autoincrement=True, primary_key=True)
    name = Column(VARCHAR(100), unique=True)
    
    position = Column(Enum(['anywhere','c-terminal','n-terminal','core']))
    
    accession = Column(VARCHAR(10))
    target = Column(VARCHAR(1))
    mono_mass_diff = Column(Float)
    avg_mass_diff = Column(Float)

    parent_id = Column(Integer(10), ForeignKey('PTM.id'))
    children = relationship("PTM")

    taxons = relationship("Taxonomy", secondary=PTM_taxon)
    keywords = relationship("PTMkeyword")

    def hasTaxon(self, taxon_ids):
        return len(set([t.node_id for t in self.taxons]) & taxon_ids) > 0

    def hasTarget(self, residue):
        return residue in set([c.target for c in self.children]) or residue == self.target
    
    def hasKeyword(self, key):
        k = key.lower()
        return k in set([kw.keyword.lower() for kw in self.keywords])

    def createKeyword(self, key):
        if not self.hasKeyword(key):
            ptmkw = PTMkeyword()
            ptmkw.keyword = key
            self.keywords.append(ptmkw)
            
class ScansitePrediction(Base):
    __tablename__ = 'phosphopep_prediction'
    id = Column(Integer(10), autoincrement=True, primary_key=True)
    source = Column(VARCHAR(40), default='scansite')
    value = Column(VARCHAR(20))
    score = Column(Float)
    phosphopep_id = Column(Integer(10), ForeignKey('phosphopep.id'))
    
    UniqueConstraint('source', 'value', 'phosphopep_id', name="UNIQUE_pepId_source_value")

class Phosphopep(Base):
    __tablename__ = 'phosphopep'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    pep_tryps = Column(VARCHAR(100))
    pep_aligned = Column(VARCHAR(15))
    
    site_pos = Column(Integer(10))
    site_type = Column(CHAR(1))
    modification_id = Column(Integer(10), ForeignKey('PTM.id'))
    
    pfam_site = Column(VARCHAR(45))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
    predictions = relationship(ScansitePrediction)
    
    def getPeptide(self):
        return self.pep_aligned
    
    def getName(self):
        return self.site_type + str(self.site_pos)

class MeasuredPeptide(Base):
    __tablename__ = 'MS'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    phosphopep = Column(VARCHAR(150))
    
    experiment = relationship("Experiment")
    protein = relationship("Protein")
    
    phosphopeps = relationship("Phosphopep", secondary=MS_phosphopep)
    data = relationship("ExperimentData")

    def save(self):
        DBSession.add(self)
        
class NoSuchModification(Exception):
    def __init__(self, pep_site, pep_type, protein_id):
        self.site = pep_site
        self.t = pep_type
        self.protein_id = protein_id
        
    def __repr__(self):
        return "No such peptide modification at site %d, residue: %s, protein: %d" % (self.site, self.t, self.protein_id)

def getMeasuredPeptidesByProtein(pid, user):
    modifications = DBSession.query(MeasuredPeptide).filter_by(protein_id=pid).all()
    return [ mod for mod in modifications if mod.experiment.checkPermissions(user) and mod.experiment.ready() ]

def getMeasuredPeptidesByExperiment(eid, user, pids = None):
    if(pids != None):
        modifications = DBSession.query(MeasuredPeptide).filter(and_(MeasuredPeptide.experiment_id==eid, MeasuredPeptide.protein_id.in_(pids))).all()
    else:
        modifications = DBSession.query(MeasuredPeptide).filter_by(experiment_id=eid).all()
    return [ mod for mod in modifications if mod.experiment.checkPermissions(user) and mod.experiment.ready() ]


def getModificationBySite(pep_site, pep_type, prot_id, mod_id=None):
    
    if mod_id==None:
        mod = DBSession.query(Phosphopep).filter_by(site_pos=pep_site, site_type=pep_type, protein_id=prot_id).first()
    else:
        mod = DBSession.query(Phosphopep).filter_by(site_pos=pep_site, site_type=pep_type, protein_id=prot_id, modification_id=mod_id).first()
    
    if mod == None: 
        raise NoSuchModification(pep_site, pep_type, prot_id)
    
    return mod


def findMatchingPTM(mod_type, residue=None, taxon_ids=None):
    mods = DBSession.query(PTM).join(PTMkeyword).filter(or_(PTM.accession==mod_type, PTM.name==mod_type, PTMkeyword.keyword==mod_type)).all()
    
    mods_exist = len(mods) > 0
    
    if residue:
        mods = [mod for mod in mods if mod.hasTarget(residue)]
    if taxon_ids:
        mods = [mod for mod in mods if mod.hasTaxon(taxon_ids)]
    
    return mods, mods_exist