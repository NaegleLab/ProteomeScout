from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column
from sqlalchemy.types import Integer, VARCHAR, Text
from sqlalchemy.orm import relationship
from sqlalchemy.sql.expression import and_

class SwissprotRecord(Base):
    __tablename__ = 'uniprot_swissprot'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    accession = Column(VARCHAR(20))
    locus = Column(VARCHAR(30))
    name = Column(VARCHAR(150))
    species = Column(VARCHAR(100))
    sequence = Column(Text)

    def __init__(self, name, accession, locus, species, sequence):
        self.name = name
        self.accession = accession
        self.locus = locus
        self.species = species
        self.sequence = sequence.upper()

def findPeptide(seq, species=None):
    filter_clause = SwissprotRecord.sequence.like('%' + seq + '%')

    if species:
        filter_clause = and_(filter_clause, SwissprotRecord.species==species)

    return DBSession.query(SwissprotRecord).filter(filter_clause).all()
