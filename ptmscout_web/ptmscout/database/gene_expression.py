from ptmscout.database import Base
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import VARCHAR, Text, Enum, Integer, Float
from sqlalchemy.orm import relationship

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
    
    probeset_id = Column(Integer(10), ForeignKey('expression_ann.id'))
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
    probeset_id = Column(Integer(10), ForeignKey('expression_ann.id'))
    
    def __init__(self, t, value, probeset_id):
        self.type = t
        self.value = value
        self.probeset_id = probeset_id

class ExpressionProbeset(Base):
    __tablename__ = 'expression_ann'
    id          = Column(Integer(10), primary_key=True, autoincrement=True)
    probeset_id = Column(VARCHAR(45), unique=True, index=True)
    genechip    = Column(Enum('gnf1h', 'gnf1m', 'HG-U133A'), default='gnf1h')
    species_id  = Column(Integer(10), ForeignKey('species.id'))
    symbol      = Column(VARCHAR(45), index=True)
    refseq      = Column(Text)
    uniprot     = Column(Text)
    aliases     = Column(Text)
    name        = Column(Text)
    
    species = relationship("Species")
    
    samples = relationship(ExpressionSample, lazy='subquery')
    accessions = relationship(ExpressionAccession)