from HTMLParser import HTMLParser
import json
from DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications, taxonomies
from paste.deploy.loadwsgi import appconfig
import os
import traceback
import re


class TableParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.isTd = False
        self.isTr = False
        self.isHeader = False
        self.table = []
        self.header = []
        self.accum = None
        
    
    def handle_starttag(self, tag, attrs):
        if tag.lower()=='tr':
            self.isTr = True
            self.table.append([])
            self.isHeader = False
        if tag.lower()=='th':
            self.isHeader = True
        if tag.lower()=='td' or tag.lower()=='th':
            self.isTd = True
            self.accum = ""
        if tag.lower()=='br' and self.isTd:
            self.accum += " "
            
        
    def handle_endtag(self, tag):
        if tag.lower()=='tr':
            self.isTr = False
            if self.isHeader:
                self.header = self.table.pop()
        if tag.lower()=='td' or tag.lower()=='th':
            self.isTd = False
            
            stripped = self.accum.strip()
            if stripped.find("\n\n") > -1:
                appender = [ s.strip() for s in stripped.split("\n\n") ]
            else:
                appender = stripped
                
            self.table[-1].append(appender)
            
            self.accum = None
        
    def handle_data(self, data):
        if self.isTd:
            self.accum += data
            

def load_ptmfile(ptmfile, aminocodes):
    ptm_dict = {}
    
    ptm = modifications.PTM()
    save = False
    target = ""
    for line in ptmfile:
        if line.strip() == "//":
            if save:
                if ptm.target != None:
                    ptm_dict[ptm.name.lower()] = ptm
                if ptm.target == None:
                    print "Error with ptm target:       " + ptm.name + " --- TG: " + target
                if ptm.position == None:
                    print "Warning ptm has no position: " + ptm.name
            
            
            
            ptm = modifications.PTM()
        try:
            param, value = line.strip().split("   ")
        except ValueError:
            continue
        
        value = value.strip()
        
        if param == "ID":
            ptm.name = value
            save = False
        if param == "PP":
            haystack = value.lower()
            if haystack.find('anywhere') > -1:
                ptm.position = 'anywhere'
            elif haystack.find('core') > -1:
                ptm.position = 'core'
            elif haystack.find('c-terminal') > -1:
                ptm.position = 'c-terminal'
            elif haystack.find('n-terminal') > -1:
                ptm.position = 'n-terminal'
                
        if param == "MM":
            ptm.mono_mass_diff = float(value)
        if param == "MA":
            ptm.avg_mass_diff = float(value)
        if param == "TR":
            m = re.match('^([A-Za-z ]*); (.*)', value)
            
            kingdom = m.group(1)
            taxons = m.group(2).split(",")
            
            for taxon in taxons:
                m = re.search('^taxId:([0-9]+)(?: \((.*)\))?', taxon.strip())
                
                tax_id = int(m.group(1))
                tax_name = m.group(2)
                tx = taxonomies.getTaxonomyById(tax_id)
                
                if tx != None:
                    ptm.taxons.append(tx)
                #else:
                #    print "No taxon id %d for species %s" % (tax_id, tax_name)
        
        if param == "TG":
            value = value.strip(".").lower()
            if value in aminocodes:
                ptm.target = aminocodes[value][1]
            else:
                target = value
        if param == "KW":
            
            for kw in value.lower().split(";"):
                kw = kw.strip().strip(".")
                ptmkw = modifications.PTMkeyword()
                ptmkw.keyword = kw
                ptm.keywords.append(ptmkw)
                
        if param == "AC":
            ptm.accession = value
        if param == "FT" and value != "CROSSLINK" and value != "CROSSLNK":
            save = True
    return ptm_dict



if __name__ == "__main__":
    parser = TableParser()
    
    with open('data/ptmtable.html', 'r') as htmlfile:
        parser.feed(htmlfile.read())
    
    with open('data/aminocodes.txt','r') as aminofile:
        amino_codes = json.loads(aminofile.read())
        
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'development.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        
        with open('data/ptmlist.txt', 'r') as ptmfile:
            ptm_dict = load_ptmfile(ptmfile, amino_codes)
        
        print "Loaded %d PTMs" % len(ptm_dict)
        
        print "Creating keyword map..."
        
        kw_map = {}
        for name in ptm_dict:
            for kw in ptm_dict[name].keywords:
                kw_lower = kw.keyword.lower()
                kws = kw_map.get(kw_lower, [])
                kws.append(ptm_dict[name])
                kw_map[kw_lower] = kws
        
        name_map = {}
        child_map = {}
        
        print "Adding keywords..."
        
        
        for row in parser.table:
            name = row[0]
            lower_name = row[0].lower()
            abbrev = row[1]
            mono_mass = float(row[2])
            avg_mass = float(row[3])
            children = row[4]
                                
            if children == '':
                children = []
            if isinstance(children, str):
                children = [children]
            
            name_map[lower_name] = abbrev
                                
            for child in children:
                child_map[child.lower()] = name
        
        for name in ptm_dict:
            ptm = ptm_dict[name]
            
            for k in ptm.keywords:
                key = k.keyword.lower()
                if key in name_map:
                    ptm.createKeyword(name_map[key])
            
            if name in name_map:
                ptm.createKeyword(name_map[name])
                
            if name in child_map:
                k = child_map[name].lower()
                ptm.createKeyword(k)
                ptm.createKeyword(name_map[k])
        
        print "Mapping children..."
        
        for row in parser.table:
            name = row[0]
            lower_name = row[0].lower()
            abbrev = row[1]
            mono_mass = float(row[2])
            avg_mass = float(row[3])
            children = row[4]
            
            if children == '':
                children = []
            if isinstance(children, str):
                children = [children]
            
            if len(children) == 0 and lower_name not in ptm_dict and lower_name not in kw_map:
                print "PTM Not Found: " + name
            if len(children) > 0:
                if lower_name not in ptm_dict:
                    ptm = modifications.PTM()
                    ptm.name = name
                    ptm.position = None
                    ptm.target = None
                    ptm.avg_mass_diff = avg_mass
                    ptm.mono_mass_diff = mono_mass
                    ptm.createKeyword(abbrev)
                    ptm_dict[lower_name] = ptm
                
                for child in children:
                    child_lower = child.lower()
                    if child_lower in ptm_dict:
                        ptm.children.append(ptm_dict[child_lower])
                    else:
                        print "Child PTM not found: " + child
        
        #first = [ ptm_dict[name] for name in ptm_dict if ptm_dict[name].parent_id == None ]
        #second = [ ptm_dict[name] for name in ptm_dict if ptm_dict[name].parent_id != None ]
        
        print "Saving session..."
        for name in ptm_dict:
            DBSession.add(ptm_dict[name])
            
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()