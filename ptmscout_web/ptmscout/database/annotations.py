from ptmscout.database import Base, DBSession
from sqlalchemy.schema import Column, ForeignKey
from sqlalchemy.types import Integer, VARCHAR, PickleType, Enum, Text
from sqlalchemy.orm import relationship


class AnnotationPermission(Base):
    __tablename__="annotation_permissions"
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    user_id = Column(Integer(10), ForeignKey('users.id'))
    annotation_set_id = Column(Integer(10), ForeignKey('annotation_sets.id'))
    access_level = Column('access_level', Enum(['view', 'owner']), default='view')

    annotation_set = relationship("AnnotationSet", backref='permissions')

    def save(self):
        DBSession.add(self)
        DBSession.flush()
class Annotation(Base):
    __tablename__ = 'MS_annotations'

    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    MS_id = Column(Integer(10), ForeignKey('MS.id'))
    type_id = Column(Integer(10), ForeignKey('annotations.id'))
    
    value = Column(VARCHAR(100), nullable=True)

class AnnotationType(Base):
    __tablename__ = 'annotations'
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    
    set_id = Column(Integer(10), ForeignKey('annotation_sets.id'))

    name = Column(VARCHAR(100), index=True)
    type = Column(Enum(['cluster', 'numeric', 'nominative']))

    order = Column(Integer(10))
    
    annotations = relationship(Annotation)
    
    def save(self):
        DBSession.add(self)
        DBSession.flush()
        
    def delete(self):
        DBSession.delete(self)

class AnnotationSet(Base):
    __tablename__ = 'annotation_sets'
    id = Column(Integer(10), primary_key=True, autoincrement=True)

    name = Column(VARCHAR(200))

    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))
    annotation_types = relationship("AnnotationType")

    def save(self):
        DBSession.add(self)
        DBSession.flush()

    def delete(self):
        DBSession.delete(self)

def getUserAnnotations(annotation_set_id, exp_id, user):
    return DBSession.query(AnnotationSet).join(AnnotationPermission).filter(AnnotationSet.experiment_id==exp_id, AnnotationPermission.user_id==user.id, AnnotationSet.id==annotation_set_id).first()

def getUserAnnotationSets(exp_id, user):
    return DBSession.query(AnnotationSet).join(AnnotationPermission).filter(AnnotationSet.experiment_id==exp_id, AnnotationPermission.user_id==user.id).all()

def getAnnotationValues(annotation_type_id):
    return DBSession.query(Annotation.value).filter_by(type_id=annotation_type_id).distinct()
    
class Subset(Base):
    __tablename__ = 'experiment_subsets'
    
    id = Column(Integer(10), primary_key=True, autoincrement=True)
    owner_id = Column(Integer(10), ForeignKey('users.id'))
    experiment_id = Column(Integer(10), ForeignKey('experiment.id'))
    annotation_set_id = Column(Integer(10), ForeignKey('annotation_sets.id'))

    name = Column(VARCHAR(100), index=True)
    
    foreground_query = Column(PickleType)
    background_query = Column(PickleType)

    share_token = Column(Text)
    
    def save(self):
        DBSession.add(self)
        DBSession.flush()
        
    def delete(self):
        DBSession.delete(self)

    def copy(self):
        c = Subset()
        c.owner_id = self.owner_id
        c.experiment_id = self.experiment_id
        c.annotation_set_id = self.annotation_set_id

        c.name = self.name

        c.foreground_query = self.foreground_query
        c.background_query = self.background_query
        c.share_token = None
        return c

def getSubsetByShareToken(share_token):
    return DBSession.query(Subset).filter_by( share_token = share_token).first()

def getSubsetByName(exp_id, subset_name, user):
    return DBSession.query(Subset).filter_by( experiment_id=exp_id, owner_id=user.id, name=subset_name ).first()

def getSubsetById(subset_id, exp_id):
    return DBSession.query(Subset).filter_by(id=subset_id, experiment_id=exp_id).first()

def getSubsetsForUser(exp_id, user):
    return DBSession.query(Subset).filter_by(experiment_id=exp_id, owner_id=user.id).all()
