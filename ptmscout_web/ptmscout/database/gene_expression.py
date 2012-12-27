from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import VARCHAR, Text, Enum, Integer, Float
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import and_

class ExpressionCollection(Base):
    __tablename__ = 'expression_collection'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    name = Column(VARCHAR(45), unique=True)

class ExpressionTissue(Base):
    __tablename__ = 'expression_tissue'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    name = Column(VARCHAR(45), unique=True)

class ExpressionSample(Base):
    __tablename__ = 'expression_samples'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    probeset_id = Column(Integer(10), ForeignKey('expression.id'))
    collection_id = Column(Integer(10), ForeignKey('expression_collection.id'))
    tissue_id = Column(Integer(10), ForeignKey('expression_tissue.id'))
    
    value = Column(Float)
    
    collection = relationship(ExpressionCollection)
    tissue = relationship(ExpressionTissue)

class ExpressionAccession(Base):
    __tablename__ = 'expression_acc'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    type = Column(Enum(['gene_symbol','refseq','uniprot', 'alias']))
    value = Column(VARCHAR(45))
    probeset_id = Column(Integer(10), ForeignKey('expression.id'))
    
    def __init__(self, t, value, probeset_id):
        self.type = t
        self.value = value
        self.probeset_id = probeset_id

class ExpressionProbeset(Base):
    __tablename__ = 'expression'
    id          = Column(Integer(10), primary_key=True, autoincrement=True)
    probeset_id = Column(VARCHAR(45), unique=True, index=True)
    genechip    = Column(Enum('gnf1h', 'gnf1m', 'HG-U133A'), default='gnf1h')
    species_id  = Column(Integer(10), ForeignKey('species.id'))
    name        = Column(Text)
    
    species = relationship("Species")
    
    samples = relationship(ExpressionSample, lazy='subquery')
    accessions = relationship(ExpressionAccession, lazy='subquery')
    
    
    
def getExpressionProbeSetsForProtein(protein_accessions, species_id):
    probesets = \
        DBSession \
            .query(ExpressionProbeset) \
            .join(ExpressionAccession) \
            .filter(
                and_(
                    ExpressionAccession.value.in_(protein_accessions),
                    ExpressionProbeset.species_id==species_id)).all()

    return probesets
    
    
