import datetime
import re

class GeneOntologyTerm():
    def __init__(self):
        self.id = None
        self.name = None
        self.namespace = None
        self.is_a_relationships = {}
        
    def get_aspect(self):
        if self.namespace == 'biological_process':
            return 'P'
        if self.namespace == 'molecular_function':
            return 'F'
        if self.namespace == 'cellular_component':
            return 'C'

    def add_is_a_relationship(self, go_id, name):
        if go_id in self.is_a_relationships:
            raise KeyError("Key: " + go_id + " already mapped")
        self.is_a_relationships[go_id] = name
        
        
class GeneOntologyFile():
    def __init__(self, filename):
        self.go_file = open(filename, 'r')
        self.last_line = "None"

        self.scan_header()
    
    def scan_header(self):
        while self.last_line != "":
            m = re.match(r"format-version: ([0-9\.]+)", self.last_line)
            if m:
                self.format_version = m.group(1)

            m = re.match(r"date: ([0-9]+):([0-9]+):([0-9]+) ([0-9]+):([0-9]+)", self.last_line)
            if m:
                year = int(m.group(3))
                month = int(m.group(2))
                day = int(m.group(1))
                hour = int(m.group(4))
                minute = int(m.group(5))
                self.date = datetime.datetime(year, month, day, hour, minute)

            self.last_line = self.go_file.readline().strip()
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


