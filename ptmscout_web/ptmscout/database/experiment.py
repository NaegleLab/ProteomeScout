from . import Base, DBSession
from sqlalchemy import Column, Integer, VARCHAR, Text

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
    journal = Column(VARCHAR(45))
    pub_date = Column(VARCHAR(20))
    
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

class NoSuchExperiment(Exception):
    def __init__(self, eid):
        self.eid = eid
    
    def __str__(self):
        return "No such experiment: %d" % (self.eid)
    
def getExperimentById(experiment_id):
    value = DBSession.query(Experiment).filter_by(id=experiment_id).first()
    if value == None:
        raise NoSuchExperiment(experiment_id)
    return value

def getAllExperiments():
    return DBSession.query(Experiment).all()

def getExperimentTree():
    experiments = getAllExperiments()
    for experiment in experiments:
        experiment.children = []
        
    experiment_dict = dict([(experiment.id, experiment) for experiment in experiments])
    experiment_tree = []
    for experiment in experiments:
        if experiment.experiment_id != 0:
            experiment_dict[experiment.experiment_id].children.append(experiment)
        else:
            experiment_tree.append(experiment)
    
    return experiment_tree
    