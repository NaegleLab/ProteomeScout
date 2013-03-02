from ptmscout.database import Base
from sqlalchemy.types import VARCHAR, Integer, Enum, Text
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.orm import relationship

class Job(Base):
    __tablename__ = 'jobs'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    status = Column(Enum(['configuration', 'in queue', 'started', 'finished', 'error']), default='configuration')
    failure_reason = Column(Text, default="")
    
    stage = Column(VARCHAR(20))
    progress = Column(Integer(10))
    max_progress = Column(Integer(10))
    
    status_url = Column(VARCHAR(250))
    resume_url = Column(VARCHAR(250))
    result_url = Column(VARCHAR(250))
    
    job_name = Column(VARCHAR(250))
    
    user_id = Column(Integer(10), ForeignKey('users.id'))
    user = relationship("User")
