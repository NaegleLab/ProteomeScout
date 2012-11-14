from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, Enum, VARCHAR
from sqlalchemy.orm import relationship

class Permission(Base):
    __tablename__="permissions"
    user_id = Column(Integer(10), ForeignKey('users.id'), primary_key=True)
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'), primary_key=True)
    access_level = Column('access_level', Enum(['view', 'owner']), default='view')
    
    def __init__(self, experiment, access_level='view'):
        self.experiment = experiment
        self.access_level = access_level

class Invitation(Base):
    __tablename__ = "invitations"
    invited_email = Column(VARCHAR(50), primary_key=True)
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'), primary_key=True)
    inviting_user_id = Column(Integer(10), ForeignKey('users.id'))

    experiment = relationship("Experiment")
    
    def __init__(self, invited_email, experiment_id, inviting_user_id):
        self.invited_email = invited_email
        self.experiment_id = experiment_id
        self.inviting_user_id = inviting_user_id
        
    def saveInvitation(self):
        DBSession.add(self)


def getInvitationsForUser(email):
    return DBSession.query(Invitation).filter_by(invited_email=email).all()