from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, Text
from ptmscout.config import settings as config
from ptmscout.database.permissions import Permission
from sqlalchemy.schema import ForeignKey, Table
from sqlalchemy.types import Float, Enum, DateTime
from sqlalchemy.sql.expression import  or_, and_
import datetime
from sqlalchemy.orm import relationship
import logging

log = logging.getLogger('ptmscout')

experiment_PTM = Table('experiment_PTM', Base.metadata,
                    Column('id', Integer(10), primary_key=True, autoincrement=True),
                    Column('experiment_id', Integer(10), ForeignKey('experiment.id')),
                    Column('PTM_id', Integer(10), ForeignKey('PTM.id')))

class PermissionException(Exception):
    def __init__(self, msg=''):
        self.msg = msg

    def __repr__(self):
        return self.msg

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
    __tablename__ = 'MS_data'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    type = Column(Enum(['data','stddev']), default='data')
    units = Column(VARCHAR(20), default='time')
    run = Column(VARCHAR(20), default='average')
    label = Column(VARCHAR(45), default='')
    
    priority = Column(Integer(10), default=0)
    value = Column(Float)
    
    MS_id = Column(Integer(10), ForeignKey('MS.id'))
    MS = relationship("MeasuredPeptide")
        
    def __format_name(self):
        return "%s:%s:%s" % (self.run, self.type, self.label)
    
    formatted_label = property(__format_name)
    
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
    experiment_id = Column(Integer(10), ForeignKey("experiment.id"), default=None)
    
    dataset = Column(Text)
    
    volume = Column(Integer(11))
    page_start = Column(VARCHAR(10))
    page_end = Column(VARCHAR(10))
    journal = Column(VARCHAR(45))
    publication_year = Column(Integer(4))
    publication_month = Column(Enum(['','january','february','march','april','may','june','july','august','september','october','november','december']))
    
    public = Column(Integer(1), default=0)
    
    job_id = Column(Integer(10))
    submitter_id = Column(Integer(10), ForeignKey('users.id'))

    modified_residues = Column(VARCHAR(40), default="")
    type = Column(Enum(['compendia','experiment','dataset'], default='experiment'))
  
    version_number = Column(Integer(11))

    def __repr__(self):
        return 'experiment:%d:%d' % (self.id, self.version_number)

    def __get_job(self):
        from ptmscout.database import jobs

        if self.job_id == None:
            return None

        if not hasattr(self, '__job_obj'):
            self.__job_obj = jobs.getJobById(self.job_id)
        
        return self.__job_obj
        
    job = property(__get_job)
    
    errors = relationship("ExperimentError", cascade="all,delete-orphan")
    conditions = relationship("ExperimentCondition", cascade="all,delete-orphan")
    permissions = relationship("Permission", backref="experiment", cascade="all,delete-orphan")
    measurements = relationship("MeasuredPeptide")
    modifications = relationship("PTM", secondary=experiment_PTM)

    parent = relationship("Experiment", backref='children', remote_side=[id])
    
    def __init__(self):
        self.date = datetime.datetime.now()
        self.version_number = 0

    def saveExperiment(self):
        DBSession.add(self)
        self.version_number += 1
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
        self.experiment_id = exp.experiment_id
        
        self.dataset = exp.dataset
        
        self.volume = exp.volume
        self.page_start = exp.page_start
        self.page_end = exp.page_end
        self.journal = exp.journal
        self.publication_year = exp.publication_year
        self.publication_month = exp.publication_month
        
        self.public = exp.public
        self.type = exp.type
        
        self.submitter_id = exp.submitter_id
        self.job_id = None
        
        self.conditions = []
        for c in exp.conditions:
            expc = ExperimentCondition()
            expc.type = c.type
            expc.value = c.value
            
            self.conditions.append(expc)
            
    def isExperiment(self):
        return self.type in set(['compendia','experiment'])
            
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
                    string += "-" + self.publication_month.capitalize()
                string += ". "
            if self.volume != "" and self.volume != None:
                string += "Vol " + str(self.volume) + ". "
            if self.page_start != None:
                string += self.page_start

                if self.page_end != None and self.page_end != self.page_start:
                    if len(self.page_start) == len(self.page_end):
                        string += "-" + self.__getPageEndFormatted()
                    else:
                        string += "-" + self.page_end

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
                    cite_string += "-" + self.publication_month.capitalize()
                cite_string += ". "
            if self.volume != None:
                cite_string += "Vol " +str(self.volume)+". "
            if self.page_start != None:
                cite_string += self.page_start

                if self.page_end != None and self.page_end != self.page_start:
                    if len(self.page_start) == len(self.page_end):
                        cite_string += "-" + self.__getPageEndFormatted()
                    else:
                        cite_string += "-" + self.page_end

                cite_string += "."
        return cite_string
   
    def hasPMID(self):
        return self.PMID != None

    def getPubMedUrl(self):
        return config.pubmedUrl % (self.PMID)

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

    def revokePermission(self, user):
        perm = None
        for p in self.permissions:
            if p.user_id == user.id:
                perm = p
        if perm != None:
            self.permissions.remove(perm)

    def __get_status(self):
        if self.job == None:
            return 'configuration'
        return self.job.status
    
    def __get_stage(self):
        if self.job == None:
            return ''
        return self.job.stage
    
    status = property(__get_status)
    loading_stage = property(__get_stage)
    
    def ready(self):
        if self.job != None:
            return self.job.status == 'finished'  
        return False
        
    def makePublic(self):
        if self.experiment_id != None:
            if self.parent.public != 1:
                raise PermissionException()
        self.public = 1
        
    def makePrivate(self):
        for child in self.children:
            if child.public == 1:
                raise PermissionException()
        self.public = 0
        
    def isOwner(self, user):
        if user == None:
            return False
        
        owner_users = [ p.user_id for p in self.permissions if p.access_level == 'owner' ]
        return user.id in owner_users
        
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

    def getClickable(self, request, new_window=False):
        url = None
        if self.type == 'experiment':
            url = request.route_url('experiment', id=self.id)
        elif self.URL and self.URL != '':
            url = self.URL

        if url:
            return '<a href="%s" %s>%s</a>' % ( url, 'target="_blank"' if new_window else '', self.name)

        return self.name

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

def getAllExperiments(current_user, filter_compendia=True, filter_datasets=True):
    exps = [ exp for exp in DBSession.query(Experiment).all() if exp.checkPermissions(current_user) and exp.ready() ]

    if filter_compendia:
        exps = [exp for exp in exps if exp.type!='compendia']

    if filter_datasets:
        exps = [exp for exp in exps if exp.type!='dataset']

    return sorted(exps, key=lambda item: item.name)

def recursive_append(exp, tree, experiments, lp=None):
    available = exp in experiments

    if available:
        tree.append(exp)
        if lp != None:
            exp.experiment_id = lp.id
        else:
            exp.experiment_id = None

        lp = exp
        experiments.remove(exp)

    for child in exp.children:
        recursive_append(child, tree, experiments, lp=lp)

def getExperimentTree(current_user):
    experiments = getAllExperiments(current_user)
    root_experiments = [ exp for exp in experiments if exp.experiment_id == None ]
    experiment_tree = []

    while len(root_experiments) > 0:
        exp = root_experiments.pop(0)
        recursive_append(exp, experiment_tree, experiments)

    return experiment_tree
    
def getValuesForField(field_name):
    results = DBSession.query(ExperimentCondition.value).filter(ExperimentCondition.type==field_name).distinct()
    return sorted([ r.value for r in results ])

def countErrorsForExperiment(exp_id):
    return DBSession.query(ExperimentError).filter_by(experiment_id=exp_id).count()

def createExperimentError(exp_id, line, accession, peptide, message):
    err = ExperimentError()
    err.experiment_id = exp_id
    err.line = line
    err.accession = accession
    err.peptide = peptide
    err.message = message
    
    DBSession.add(err)

def errorsForAccession(exp_id, accession):
    return DBSession.query(ExperimentError).filter_by(experiment_id=exp_id, accession=accession).all()

def searchExperiments(text_search=None, conditions={}, user=None, page=None):
    q = DBSession.query(Experiment.id).outerjoin(Experiment.conditions).outerjoin(Experiment.permissions)

    filter_clauses = [ Experiment.type != 'dataset' ]
    if user:
        filter_clauses.append( or_( Experiment.public==1, Permission.user_id==user.id ) )
    else:
        filter_clauses.append( Experiment.public == 1 )

    if text_search:
        text_search = '%' + str(text_search) + '%'
        filter_clauses.append( or_( Experiment.name.like(text_search), Experiment.description.like(text_search) ) )

    for cond_type in conditions:
        for cond_value in conditions[cond_type]:
            filter_clauses.append( and_( ExperimentCondition.type==cond_type, ExperimentCondition.value==cond_value ) )

    if len(filter_clauses) > 0:
        filter_clause = and_(*filter_clauses)
    else:
        filter_clause = filter_clauses[0]

    q = q.filter( filter_clause ).distinct().order_by( Experiment.name )

    if page==None:
        sq = q.subquery()
    else:
        limit, offset = page
        sq = q.limit(limit).offset(offset).subquery()

    fq = DBSession.query(Experiment).join(sq, Experiment.id == sq.c.id)
    return q.count(), fq.all()

def countExperiments():
    exps = [ exp for exp in DBSession.query(Experiment).all() if exp.ready() ]

    compendia = len([ exp for exp in exps if exp.type == 'compendia' ])
    experiments = len([ exp for exp in exps if exp.type == 'experiment' ])
    datasets = len([ exp for exp in exps if exp.type == 'dataset' ])

    return compendia, experiments, datasets

