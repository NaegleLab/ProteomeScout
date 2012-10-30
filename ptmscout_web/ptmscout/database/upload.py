from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR, Enum, Text, DateTime
from sqlalchemy.orm import relationship
import datetime

class SessionColumn(Base):
    __tablename__='session_columns'
    id=Column(Integer(10), primary_key=True, autoincrement=True)
    session_id=Column(Integer(10), ForeignKey('sessions.id'))
    
    column_values=['hidden','data','stddev','accession','peptide','species','modification','run', 'none']
    
    type=Column(Enum(column_values), default='none')
    label=Column(VARCHAR(45), default='')
    column_number=Column(Integer)

class Session(Base):
    __tablename__='sessions'
    id=Column(Integer(10), primary_key=True, autoincrement=True)
    user_id=Column(Integer(10), ForeignKey('users.id'))
    
    data_file=Column(VARCHAR(100))
    load_type=Column(Enum(['new','reload','append','extension']))
    parent_experiment=Column(Integer(10), ForeignKey('experiment.id'))
    change_description=Column(Text)
    units=Column(VARCHAR(20), default='')
    
    stage=Column(Enum(['config','metadata','confirm','complete']))
    experiment_id=Column(Integer(10), ForeignKey('experiment.id'))
    date=Column(DateTime)
    
    columns = relationship("SessionColumn")
        
    def __init__(self):
        self.date = datetime.datetime.now()
    
    def save(self):
        DBSession.add(self)
        DBSession.flush()
        
    def getAncestor(self):
        ancestors = DBSession.query(Session).filter(Session.experiment_id==self.parent_experiment).all()
        ancestors.sort(key=lambda session: session.date, reverse=True)
        return ancestors[0]

class NoSuchSession(Exception):
    pass


class SessionAccessForbidden(Exception):
    pass


def getSessionById(session_id, current_user):
    session = DBSession.query(Session).filter_by(id=session_id).first()
    
    if session == None:
        raise NoSuchSession()
    
    if current_user == None or session.user_id != current_user.id:
        raise SessionAccessForbidden()
         
    
    return session