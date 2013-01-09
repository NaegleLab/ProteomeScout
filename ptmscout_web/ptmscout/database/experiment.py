from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, Text
from ptmscout.config import settings as config
from ptmscout.database.permissions import Permission
from sqlalchemy.schema import ForeignKey
from sqlalchemy.types import Float, Enum, DateTime
from sqlalchemy.sql.expression import null
import datetime
from sqlalchemy.orm import relationship, deferred
import logging
import pickle

log = logging.getLogger('ptmscout')

class ExperimentError(Base):
    __tablename__ = 'experiment_error'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    experiment_id = Column(Integer(10), ForeignKey("experiment.id"))
    line = Column(Integer(10))
    
    accession = Column(VARCHAR(45))
    peptide = Column(VARCHAR(150))
    
    message = Column(Text)
    
    
class ExperimentCondition(Base):
    __tablename__ = 'experiment_condition'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    type = Column(Enum(['cell','tissue','drug','stimulus','environment']))
    value = Column(VARCHAR(100))
    
    experiment_id = Column(Integer(10), ForeignKey("experiment.id"))

class ExperimentData(Base):
    __tablename__ = 'data'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    type = Column(Enum(['data','stddev']), default='data')
    units = Column(VARCHAR(20), default='time')
    run = Column(VARCHAR(20), default='average')
    label = Column(VARCHAR(45), default='')
    
    priority = Column(Integer(10), default=0)
    value = Column(Float, default=null)
    
    MS_id = Column(Integer(10), ForeignKey('MS.id'))
    MS = relationship("MeasuredPeptide")
    
    def save(self):
        DBSession.add(self)

class Experiment(Base):
    __tablename__ = 'experiment'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    name = Column(Text)
    author = Column(Text)
    date = Column(DateTime)
    
    description = Column(Text)
    contact = Column(VARCHAR(45))
    
    PMID = Column(Integer(10))
    URL = Column(Text)
    published = Column(Integer(1))
    
    ambiguity = Column(Integer(1))
    export = Column(Integer(1))
    experiment_id = Column(Integer(10), default=None)
    
    dataset = Column(Text)
    
    volume = Column(Integer(11))
    page_start = Column(VARCHAR(10))
    page_end = Column(VARCHAR(10))
    journal = Column(VARCHAR(45))
    publication_year = Column(Integer(4))
    publication_month = Column(Enum(['','january','february','march','april','may','june','july','august','september','october','november','december']))
    
    public = Column(Integer(1), default=0)
    
    status = Column(Enum('configuration', 'in queue','loading','loaded', 'error'), default='configuration')
    submitter_id = Column(Integer(10), ForeignKey('users.id'))

    loading_stage = Column(Enum('query', 'proteins', 'GO terms', 'peptides', 'scansite'), default='query')
    progress = Column(Integer(10), default=0)
    max_progress = Column(Integer(10), default=0)

    last_stage_result = deferred(Column(Text, default="N."))
    failure_reason = Column(Text, default="")
    
    errors = relationship("ExperimentError", cascade="all,delete-orphan")
    conditions = relationship("ExperimentCondition", cascade="all")
    permissions = relationship("Permission", backref="experiment", cascade="all,delete-orphan")
    
    def __init__(self):
        self.date = datetime.datetime.now()

    def store_last_result(self, result):
        self.last_stage_result = pickle.dumps(result)

    def get_last_result(self):
        try:
            return pickle.loads(self.last_stage_result)
        except:
            return None

    def saveExperiment(self):
        DBSession.add(self)
        DBSession.flush()
        
    def delete(self):
        DBSession.delete(self)
        
    def copyData(self, exp):
        self.name = exp.name
        self.author = exp.author
        self.date = exp.date
        
        self.description = exp.description
        self.contact = exp.contact
        
        self.PMID = exp.PMID
        self.URL = exp.URL
        self.published = exp.published
        
        self.ambiguity = exp.ambiguity
        self.export = exp.export
        self.experiment_id = exp.experiment_id
        
        self.dataset = exp.dataset
        
        self.volume = exp.volume
        self.page_start = exp.page_start
        self.page_end = exp.page_end
        self.journal = exp.journal
        self.publication_year = exp.publication_year
        self.publication_month = exp.publication_month
        
        self.public = exp.public
        
        self.status = exp.status
        self.submitter_id = exp.submitter_id
        
        self.loading_stage = 'query'
        self.progress = 0
        self.max_progress = 0
        self.last_stage_result = ''
        self.failure_reason = ''
        
        self.conditions = []
        for c in exp.conditions:
            expc = ExperimentCondition()
            expc.type = c.type
            expc.value = c.value
            
            self.conditions.append(expc)
            
    def clearErrors(self):
        self.errors = []

    def getCitationString(self):
        string = ""
        if self.published == 1:
            if self.journal != "" and self.journal != None:
                string += self.journal + ". "
            if self.publication_year != None:
                string += str(self.publication_year)
                if self.publication_month != None:
                    string += "-" + self.publication_month
                string += ". "
            if self.volume != "" and self.volume != None:
                string += "Vol " + str(self.volume) + ". "
            if self.page_start != None:
                string += self.page_start
                if self.page_end != None:
                    string += "-" + self.__getPageEndFormatted()
                string += "."
        return string
    
    def __getPageEndFormatted(self):
        i = 0
        while(self.page_start[i] == self.page_end[i]):
            i+=1
        return self.page_end[i:]
    
    def getLongCitationString(self):
        cite_string = ""
        if self.published == 1:
            if self.author !="":
                cite_string += self.author + ". "
            if self.journal != "":
                cite_string += "<b>"+self.journal + "</b>. "
            if self.publication_year != None:
                cite_string += str(self.publication_year)
                if self.publication_month != None:
                    cite_string += "-" + self.publication_month
                cite_string += ". "
            if self.volume != None:
                cite_string += "Vol " +str(self.volume)+". "
            if self.page_start != None:
                cite_string += self.page_start
                if self.page_end != None:
                    cite_string += "-" + self.__getPageEndFormatted()
                cite_string += "."
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
        found = False
        for p in self.permissions:
            if p.user_id == user.id:
                p.access_level = level
                found = True
        
        if not found:
            p = Permission(self, level)
            p.user = user
            p.access_level = level
            self.permissions.append(p)
    
    def ready(self):
        return self.status == 'loaded'
        
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
    
    def hasExperimentCondition(self, ctype, value):
        for cond in self.conditions:
            if cond.type == ctype and cond.value == value:
                return True    
        return False
    
    def addExperimentCondition(self, ctype, value):
        if self.hasExperimentCondition(ctype, value):
            return 
        
        condition = ExperimentCondition()
        condition.type = ctype
        condition.value = value
        self.conditions.append(condition)
        
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

class ExperimentNotAvailable(Exception):
    def __init__(self, eid):
        self.eid = eid
    
    def __str__(self):
        return "Experiment %d is still being processed for upload" % (self.eid)
    

def getExperimentById(experiment_id, user=None, check_ready=True, secure=True):
    value = DBSession.query(Experiment).filter_by(id=experiment_id).first()
    if value == None:
        raise NoSuchExperiment(experiment_id)

    if secure and not value.checkPermissions(user):
        raise ExperimentAccessForbidden(experiment_id)
    
    if check_ready and not value.ready():
        raise ExperimentNotAvailable(experiment_id)
    
    return value

def getAllExperiments(current_user):
    return [ exp for exp in  DBSession.query(Experiment).all() if exp.checkPermissions(current_user) and exp.ready() and exp.export == 1 ]

def getExperimentTree(current_user):
    experiments = getAllExperiments(current_user)
    for experiment in experiments:
        experiment.children = []
        
    experiment_dict = dict([(experiment.id, experiment) for experiment in experiments])
    experiment_tree = []
    for experiment in experiments:
        if experiment.experiment_id != None:
            if experiment.experiment_id in experiment_dict:
                experiment_dict[experiment.experiment_id].children.append(experiment)
        else:
            experiment_tree.append(experiment)
    
    return experiment_tree
    
def getValuesForField(field_name):
    results = DBSession.query(ExperimentCondition).filter(ExperimentCondition.type==field_name).all()
    
    return sorted([ r.value for r in results ])

def createExperimentError(exp_id, line, accession, peptide, message):
    log.debug("Error: (%s, %s) %s", accession, peptide, message)
    err = ExperimentError()
    err.experiment_id = exp_id
    err.line = line
    err.accession = accession
    err.peptide = peptide
    err.message = message
    
    DBSession.add(err)
    
