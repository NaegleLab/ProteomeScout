from ptmscout.database import Base
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR, DateTime
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.mysql import INTEGER as Integer 
import datetime
import re

class Mutation(Base):
    __tablename__= 'protein_mutations'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    protein_id = Column(Integer(unsigned=True, width=10), ForeignKey('protein.id'))
    acc_id = Column(VARCHAR(30))
    mutationType = Column(VARCHAR(25))
    location = Column(Integer(5))
    original = Column(VARCHAR(3))
    mutant = Column(VARCHAR(3))
    date = Column(DateTime)
    annotation = Column(VARCHAR(256))
    clinical = Column(VARCHAR(100))

    protein = relationship("Protein")

    def __init__(self, mtype, location, original, mutant, acc_ID, annotation, PID):
        self.protein_id = PID
        self.acc_id = acc_ID
        self.mutationType = mtype
        self.location = location
        self.original = original
        self.mutant = mutant
        self.date = datetime.datetime.now()
        self.annotation = annotation

    def __repr__(self):
        return "%s%d%s" % (self.original, self.location, self.mutant)

    def equals(self, other):
        type_eq = self.mutationType == other.mutationType
        loc_eq = self.location == other.location
        mut_eq = self.mutant == other.mutant

        return type_eq and loc_eq and mut_eq

    def consistent(self, prot_seq):
        start = self.location-1
        return prot_seq[start:start+len(self.original)] == self.original


    def getDBSNPId(self):
        m = re.search('rs[0-9]+', self.annotation)
        if m == None:
            return None
        return m.group(0)
