from pyramid import testing
from paste.deploy.loadwsgi import appconfig
from sqlalchemy import engine_from_config
from ptmscout.database import DBSession, Base  # base declarative object
import os
import traceback
from ptmscout.database.protein import GeneOntology

class DatabaseInitialization():
    @classmethod
    def setUpClass(cls):
        cls.engine = engine_from_config(settings, prefix='sqlalchemy.')

    def setUp(self):
        self.connection = self.engine.connect()

        # begin a non-ORM transaction
        self.trans = self.connection.begin()

        # bind an individual Session to the connection
        DBSession.configure(bind=self.connection)
        self.session = DBSession
        Base.session = self.session

    def rollback(self):
        testing.tearDown()
        self.trans.rollback()
        self.session.close()

    def tearDown(self):
        testing.tearDown()
        self.trans.commit()
        self.session.close()
        
class GeneOntologyTerm():
    def __init__(self):
        self.id = None
        self.name = None
        self.namespace = None
        self.is_a_relationships = {}
        
    def add_is_a_relationship(self, go_id, name):
        if go_id in self.is_a_relationships:
            raise KeyError("Key: " + go_id + " already mapped")
        self.is_a_relationships[go_id] = name
        
        
class GeneOntologyFile():
    def __init__(self, filename):
        self.go_file = open(filename, 'r')
        self.last_line = "None"
    
    def __to_line(self, match):
        while(self.last_line != "" and self.last_line.strip() != match):
            self.last_line = self.go_file.readline()
            
        return self.last_line.strip() == match
            
    
    def next_term(self):
        found = self.__to_line("[Term]")
        
        ID_TAG = "id:"
        NAME_TAG = "name:"
        NAMESPACE_TAG = "namespace:"
        IS_A_TAG = "is_a:"

        if not found:
            return None
        
        term = None
        
        while(self.last_line.strip() != ""):
            self.last_line = self.last_line.strip()
            
            if self.last_line == "[Term]":
                term = GeneOntologyTerm()
                
            if self.last_line.find(ID_TAG) == 0:
                term.id = self.last_line[len(ID_TAG):].strip()
                
            if self.last_line.find(NAME_TAG) == 0:
                term.name = self.last_line[len(NAME_TAG):].strip()
                
            if self.last_line.find(NAMESPACE_TAG) == 0:
                term.namespace = self.last_line[len(NAMESPACE_TAG):].strip()
                
            if self.last_line.find(IS_A_TAG) == 0:
                rel = self.last_line[len(IS_A_TAG):].strip().split("!")
                
                term.add_is_a_relationship(rel[0].strip(), rel[1].strip())

            self.last_line = self.go_file.readline()
        
        return term
        
    
    def close(self):
        self.go_file.close()

if __name__ == '__main__':
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'development.ini'))
            
        DatabaseInitialization.setUpClass()
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        
        GO_file_location = "/home/ptmscoutdev/Desktop/gene_ontology.1_2.obo"
        GO_file = GeneOntologyFile(GO_file_location)
        
        missing_terms = set()
        
        print "Loading all terms from DB..."
        GO_terms = DBSession.query(GeneOntology).all()
        
        print "mapping terms..."
        term_map = {}
        for entry in GO_terms:
            term_map[entry.GO] = entry
        
        i=0
        term = GO_file.next_term()
        print "inserting relationships..."
        while term != None:
            if(term.id not in term_map):
                missing_terms.add(term.id)
            else:
                db_entry = term_map[term.id]
                
                for pid, pterm in term.is_a_relationships.items():
                    if(pid not in term_map):
                        missing_terms.add(pid)
                    else:
                        parent_entry = DBSession.query(GeneOntology).filter_by(GO=pid).first()
                        parent_entry.children.append(db_entry)
                    
            i+=1
            if(i % 1000 == 0):
                print ".",
            term = GO_file.next_term()
        
        print "Missing terms:", len(missing_terms)
        ofile = open('missing_terms', 'w')
        for term in missing_terms:
            ofile.write(term+"\n")
        ofile.close()
        
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()