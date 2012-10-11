from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, Text
from ptmscout.config import settings as config
from ptmscout.database.permissions import Permission
from sqlalchemy.dialects.mysql.base import TINYINT
from sqlalchemy.schema import ForeignKey
from sqlalchemy.types import Float, Enum
from sqlalchemy.sql.expression import null

class ExperimentData(Base):
    __tablename__ = 'data'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    type = Column(VARCHAR(20), default='time')
    run = Column(VARCHAR(20), default='average')
    label = Column(VARCHAR(45), default='')
    
    priority = Column(Integer(10), default=0)
    value = Column(Float, default=null)
    
    NA = Column(TINYINT(1), default=0)
    MS_id = Column(Integer(10), ForeignKey('MS.id'))


class Experiment(Base):
    __tablename__ = 'experiment'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    name = Column(Text)
    author = Column(Text)
    date = Column(VARCHAR(45))
    
    description = Column(Text)
    contact = Column(VARCHAR(45))
    
    PMID = Column(Integer(10))
    URL = Column(Text)
    published = Column(Integer(1))
    
    errorLog = Column(VARCHAR(30))
    ambiguity = Column(Integer(1))
    export = Column(Integer(1))
    experiment_id = Column(Integer(10))
    
    dataset = Column(Text)
    submitter = Column(VARCHAR(50))
    
    primaryModification = Column(VARCHAR(15))
    
    volume = Column(Integer(11))
    pages = Column(VARCHAR(20))
    page_start = Column(VARCHAR(10))
    page_end = Column(VARCHAR(10))
    journal = Column(VARCHAR(45))
    pub_date = Column(VARCHAR(20))
    publication_year = Column(Integer(4))
    publication_month = Column(Enum(['','january','february','march','april','may','june','july','august','september','october','november','december']))
    
    public = Column(Integer(1),default=0)
    
    def __init__(self):
        pass

    def saveExperiment(self):
        DBSession.add(self)
        DBSession.flush()
        
    def getCitationString(self):
        string = ""
        if self.journal != "" and self.journal != None:
            string += self.journal + ". "
        if self.pub_date != "" and self.pub_date != None:
            string += self.pub_date + ". "
        if self.volume != "" and self.volume != None:
            string += "Vol " + str(self.volume) + ". "
        if self.pages != "" and self.pages != None:
            string += self.pages + "."
        return string
    
    def getLongCitationString(self):
        cite_string = ""
        if self.author !="":
            cite_string += self.author + ". "
        if self.journal != "":
            cite_string += "<b>"+self.journal + "</b>. "
        if self.pub_date != "":
            cite_string += self.pub_date + ". "
        if self.volume != None:
            cite_string += "Vol " +str(self.volume)+". "
        if self.pages != "":
            cite_string += self.pages + "."
        return cite_string
    
    def getUrl(self):
        url = self.URL
        if url == "NA":
            url = None
        
        if url != None:
            return url
        elif self.PMID != None:
            return config.pubmedUrl % (self.PMID)
        
        return None
    
    def grantPermission(self, user, level):
        p = Permission(self, level)
        p.user = user
        
        self.permissions.append(p)
        
    def makePublic(self):
        self.public = 1
        
    def makePrivate(self):
        self.public = 0
        
    def checkPermissions(self, user):
        if self.public == 1:
            return True
        
        if user == None:
            return False
        
        allowed_users = [p.user_id for p in self.permissions]
        return user.id in allowed_users
        
class NoSuchExperiment(Exception):
    def __init__(self, eid):
        self.eid = eid
    
    def __str__(self):
        return "No such experiment: %d" % (self.eid)

class ExperimentAccessForbidden(Exception):
    def __init__(self, eid):
        self.eid = eid
    
    def __str__(self):
        return "Current user does not have access privileges to experiment %d" % (self.eid)


def getExperimentById(experiment_id, current_user):
    value = DBSession.query(Experiment).filter_by(id=experiment_id).first()
    if value == None:
        raise NoSuchExperiment(experiment_id)

    if not value.checkPermissions(current_user):
        raise ExperimentAccessForbidden(experiment_id)
    
    return value

def getAllExperiments(current_user):
    return [ exp for exp in  DBSession.query(Experiment).all() if exp.checkPermissions(current_user) ]

def getExperimentTree(current_user):
    experiments = getAllExperiments(current_user)
    for experiment in experiments:
        experiment.children = []
        
    experiment_dict = dict([(experiment.id, experiment) for experiment in experiments])
    experiment_tree = []
    for experiment in experiments:
        if experiment.experiment_id != 0:
            if experiment.experiment_id in experiment_dict:
                experiment_dict[experiment.experiment_id].children.append(experiment)
        else:
            experiment_tree.append(experiment)
    
    return experiment_tree
    