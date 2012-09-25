from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey, Table
from sqlalchemy.types import Integer, VARCHAR, CHAR
from sqlalchemy.orm import relationship

MS_phosphopep = Table('MS_phosphopep', Base.metadata,
                      Column('id', Integer(10), primary_key=True),
                      Column('MS_id', Integer(10), ForeignKey('MS.id')),
                      Column('phosphopep_id', Integer(10), ForeignKey('phosphopep.id')))


class Phosphopep(Base):
    __tablename__ = 'phosphopep'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    pep_tryps = Column(VARCHAR(100))
    pep_aligned = Column(VARCHAR(15))
    site_pos = Column(Integer(10))
    site_type = Column(CHAR(1))
    pfam_site = Column(VARCHAR(45))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    
    def getPeptide(self):
        return self.pep_aligned
    
    def getName(self):
        return self.site_type + str(self.site_pos)

class Modification(Base):
    __tablename__ = 'MS'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))
    protein_id = Column(Integer(10), ForeignKey('protein.id'))
    phosphopep = Column(VARCHAR(150))
    
    experiment = relationship("Experiment")
    phosphopeps = relationship("Phosphopep", secondary=MS_phosphopep)
    data = relationship("ExperimentData")


def getModificationsByProtein(pid, user):
    modifications = DBSession.query(Modification).filter_by(protein_id=pid).all()
    return [ mod for mod in modifications if mod.experiment.checkPermissions(user) ]
    