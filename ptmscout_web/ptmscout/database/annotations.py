from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR, PickleType


#class ClusterLabel(Base):
#    __tablename__ = 'ms_cluster'
#
#    id = Column(Integer(10), primary_key=True, autoincrement=True)
#    cluster_id = Column(Integer(10), ForeignKey('clusters.id'))
#    ms_id = Column(Integer(10), ForeignKey('ms.id'))
#    
#    label = Column(VARCHAR(10), index=True)
#
#class Cluster(Base):
#    __tablename__ = 'cluster'
#    
#    id = Column(Integer(10), primary_key=True, autoincrement=True)
#    owner_id = Column(Integer(10), ForeignKey('user.id'))
#    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))
#
#    name = Column(VARCHAR(100))

class Subset(Base):
    __tablename__ = 'experiment_subsets'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    owner_id = Column(Integer(10), ForeignKey('users.id'))
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))

    name = Column(VARCHAR(100), index=True)
    
    foreground_query = Column(PickleType)
    background_query = Column(PickleType)
    
    def save(self):
        DBSession.add(self)
        DBSession.flush()
        
    def delete(self):
        DBSession.delete(self)

def getSubsetByName(exp_id, subset_name, user):
    return DBSession.query(Subset).filter_by( experiment_id=exp_id, owner_id=user.id, name=subset_name ).first()

def getSubsetById(subset_id, exp_id):
    return DBSession.query(Subset).filter_by(id=subset_id, experiment_id=exp_id).first()

def getSubsetsForUser(exp_id, user):
    return DBSession.query(Subset).filter_by(experiment_id=exp_id, owner_id=user.id).all()