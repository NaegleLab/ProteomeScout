from ptmscout.database import DBSession, Base
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, Enum
from sqlalchemy.orm import relationship

class Permission(Base):
    __tablename__="permissions"
    user_id = Column(Integer(10), ForeignKey('users.id'), primary_key=True)
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'), primary_key=True)
    access_level = Column('access_level', Enum(['view', 'owner']), default='view')
    
    experiment = relationship("Experiment", backref="permissions")

    def __init__(self, experiment, access_level='view'):
        self.experiment = experiment
        self.access_level = access_level
