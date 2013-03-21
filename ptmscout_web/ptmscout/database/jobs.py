from ptmscout.database import Base, DBSession
from sqlalchemy.types import VARCHAR, Integer, Enum, Text, DateTime
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.orm import relationship
import datetime
from ptmscout.config import settings

class Job(Base):
    __tablename__ = 'jobs'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    status = Column(Enum(['configuration', 'in queue', 'running', 'finished', 'error']), default='configuration')
    failure_reason = Column(Text, default="")
    
    stage = Column(VARCHAR(20))
    progress = Column(Integer(10))
    max_progress = Column(Integer(10))
    
    status_url = Column(VARCHAR(250), nullable=True)
    resume_url = Column(VARCHAR(250), nullable=True)
    result_url = Column(VARCHAR(250), nullable=True)
    
    name = Column(VARCHAR(250))
    type = Column(Enum(['load_experiment','load_annotations','load_dataset']))
    
    user_id = Column(Integer(10), ForeignKey('users.id'))
    user = relationship("User", backref='jobs')
    
    created = Column(DateTime)
    restarted = Column(DateTime, nullable=True)
    finished = Column(DateTime, nullable=True)

    def __init__(self):
        self.created = datetime.datetime.now()
        
        self.max_progress = 0
        self.progress = 0
        self.stage = 'initializing'
        self.status = 'configuration'
        
    def started(self):
        if self.restarted == None:
            return self.created
        return self.restarted
    
    def finished_time(self):
        if self.finished == None:
            return "-"
        return self.finished
    
    def restart(self):
        self.restarted = datetime.datetime.now()
    
    def fail(self, stack_trace):
        self.failure_reason = stack_trace
        self.status = 'error'
        self.finished = datetime.datetime.now()

    def is_active(self):
        return self.status not in ['error', 'finished']
    
    def is_old(self):
        if self.finished == None:
            return False
        delta = datetime.datetime.now() - self.finished
        return delta.total_seconds() > settings.JOB_AGE_LIMIT

    def finish(self):
        self.failure_reason = ''
        self.status = 'finished'
        self.finished = datetime.datetime.now()
        
    def save(self):
        DBSession.add(self)
        DBSession.flush()
        
        
def getJobById(jid):
    return DBSession.query(Job).filter_by(id=jid).first()
