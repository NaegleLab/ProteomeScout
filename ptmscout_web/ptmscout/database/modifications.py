from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, Table, UniqueConstraint
from sqlalchemy.types import Integer, VARCHAR, CHAR, Float, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import and_, or_

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
    parent = relationship("PTM", backref="children", remote_side='PTM.id')

    taxons = relationship("Taxonomy", secondary=PTM_taxon)
    keywords = relationship("PTMkeyword")

    def hasTaxon(self, search_taxons):
        return len(set([t.name.lower() for t in self.taxons]) & search_taxons) > 0

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
    __tablename__ = 'peptide_predictions'
    id = Column(Integer(10), autoincrement=True, primary_key=True)
    source = Column(VARCHAR(40), default='scansite')
    value = Column(VARCHAR(20))
    score = Column(Float)
    peptide_id = Column(Integer(10), ForeignKey('peptide.id'))
    
    UniqueConstraint('source', 'value', 'peptide_id', name="UNIQUE_pepId_source_value")

class Peptide(Base):
    __tablename__ = 'peptide'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    pep_aligned = Column(VARCHAR(15))
    
    site_pos = Column(Integer(10))
    site_type = Column(CHAR(1))
    
    protein_domain_id = Column(Integer(10), ForeignKey('protein_domain.id'))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
    protein_domain = relationship("ProteinDomain")
    predictions = relationship(ScansitePrediction)
    
    def getPeptide(self):
        return self.pep_aligned
    
    def getName(self):
        return self.site_type + str(self.site_pos)

    def save(self):
        DBSession.add(self)
        DBSession.flush()

class PeptideModification(Base):
    __tablename__ = 'MS_modifications'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    MS_id = Column(Integer(10), ForeignKey('MS.id'))
    peptide_id = Column(Integer(10), ForeignKey('peptide.id'))
    modification_id = Column(Integer(10), ForeignKey('PTM.id'))
    
    peptide = relationship("Peptide")
    modification = relationship("PTM")
    

class MeasuredPeptide(Base):
    __tablename__ = 'MS'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    peptide = Column(VARCHAR(150))
    
    experiment = relationship("Experiment")
    protein = relationship("Protein")
    
    peptides = relationship("PeptideModification")
    data = relationship("ExperimentData")

    def save(self):
        DBSession.add(self)
        
    def addPeptideModification(self, peptide, ptm):
        pmod = PeptideModification()
        
        pmod.peptide = peptide
        pmod.modification = ptm
        self.peptides.append(pmod)

    def hasPeptideModification(self, peptide, ptm):
        for modpep in self.peptides:
            if modpep.peptide_id == peptide.id and modpep.modification_id == ptm.id:
                return True
        return False

    def getDataElement(self, run_name, tp, x):
        for d in self.data:
            if d.run == run_name and d.type == tp and d.label == x:
                return d

class NoSuchPeptide(Exception):
    def __init__(self, pep_site, pep_type, protein_id):
        self.site = pep_site
        self.t = pep_type
        self.protein_id = protein_id
        
    def __repr__(self):
        return "No such peptide modification at site %d, residue: %s, protein: %d" % (self.site, self.t, self.protein_id)

def getMeasuredPeptide(exp_id, pep_seq, protein_id):
    return DBSession.query(MeasuredPeptide).filter_by(experiment_id=exp_id, peptide=pep_seq, protein_id=protein_id).first()


def getMeasuredPeptidesByProtein(pid, user):
    modifications = DBSession.query(MeasuredPeptide).filter_by(protein_id=pid).all()
    return [ mod for mod in modifications if mod.experiment.checkPermissions(user) and mod.experiment.ready() ]

def getMeasuredPeptidesByExperiment(eid, user=None, pids = None, secure=True, check_ready=True):
    if(pids != None):
        modifications = DBSession.query(MeasuredPeptide).filter(and_(MeasuredPeptide.experiment_id==eid, MeasuredPeptide.protein_id.in_(pids))).all()
    else:
        modifications = DBSession.query(MeasuredPeptide).filter_by(experiment_id=eid).all()
    return [ mod for mod in modifications if (not secure or mod.experiment.checkPermissions(user)) and (not check_ready or mod.experiment.ready()) ]

def countMeasuredPeptidesForExperiment(eid):
    return DBSession.query(MeasuredPeptide).filter_by(experiment_id=eid).count()

def getPeptideBySite(pep_site, pep_type, prot_id):
    mod = DBSession.query(Peptide).filter_by(site_pos=pep_site, site_type=pep_type, protein_id=prot_id).first()
    
    if mod == None: 
        raise NoSuchPeptide(pep_site, pep_type, prot_id)
    
    return mod


def findMatchingPTM(mod_type, residue=None, taxons=None):
    mods = DBSession.query(PTM).outerjoin(PTMkeyword).filter(or_(PTM.accession==mod_type, PTM.name==mod_type, PTMkeyword.keyword==mod_type)).all()
    
    mods_exist = len(mods) > 0
    
    if residue:
        mods = [mod for mod in mods if mod.hasTarget(residue)]
    if taxons:
        mods = [mod for mod in mods if mod.hasTaxon(taxons) or len(mod.taxons) == 0]
    
    return mods, mods_exist

def getPeptideById(pep_id):
    return DBSession.query(Peptide).filter(Peptide.id==pep_id).first()

def deleteExperimentData(exp_id):
    DBSession.query(MeasuredPeptide).filter_by(experiment_id=exp_id).delete()
    DBSession.flush()
