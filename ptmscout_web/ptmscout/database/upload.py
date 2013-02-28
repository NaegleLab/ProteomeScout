from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR, Enum, Text, DateTime
from sqlalchemy.orm import relationship
import datetime

class SessionColumn(Base):
    __tablename__='session_columns'
    id=Column(Integer(10), primary_key=True, autoincrement=True)
    session_id=Column(Integer(10), ForeignKey('session.id'))
    
    column_values=['hidden','data','stddev','accession','peptide','species','modification','run', 'none']
    
    type=Column(Enum(column_values), default='none')
    label=Column(VARCHAR(45), default='')
    column_number=Column(Integer)

class Session(Base):
    __tablename__='session'
    id=Column(Integer(10), primary_key=True, autoincrement=True)
    user_id=Column(Integer(10), ForeignKey('users.id'))
    
    data_file=Column(VARCHAR(100))
    load_type=Column(Enum(['new','reload','append','extension','annotations']))
    parent_experiment=Column(Integer(10), ForeignKey('experiment.id'))
    change_description=Column(Text)
    units=Column(VARCHAR(20), default='')
    
    stage=Column(Enum(['config','metadata','confirm','condition', 'complete']))
    experiment_id=Column(Integer(10), ForeignKey('experiment.id'))
    date=Column(DateTime)
    
    columns = relationship("SessionColumn", cascade="all,delete-orphan", lazy="joined")
        
    def __init__(self):
        self.date = datetime.datetime.now()
    
    def save(self):
        DBSession.add(self)
        DBSession.flush()

    def delete(self):
        DBSession.delete(self)
        
    def getAncestor(self):
        ancestors = DBSession.query(Session).filter(Session.experiment_id==self.parent_experiment).all()
        ancestors.sort(key=lambda session: session.date, reverse=True)
        return ancestors[0]

    def getColumns(self, tp):
        found_columns = []
        for col in self.columns:
            if col.type == tp:
                found_columns.append(col)
        
        return found_columns
            
class NoSuchSession(Exception):
    pass


class SessionAccessForbidden(Exception):
    pass

def getMostRecentSession(exp_id):
    ancestors = DBSession.query(Session).filter(Session.experiment_id==exp_id).all()
    ancestors.sort(key=lambda session: session.date, reverse=True)
    return ancestors[0]

def getSessionById(session_id, user=None, secure=True):
    session = DBSession.query(Session).filter_by(id=session_id).first()
    
    if session == None:
        raise NoSuchSession()
    
    if secure and (user == None or session.user_id != user.id):
        raise SessionAccessForbidden()
         
    
    return session