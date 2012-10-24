from ptmscout.database import DBSession, taxonomies
from DB_init import DatabaseInitialization
import traceback
from paste.deploy.loadwsgi import appconfig
import os
import re

def get_similar_taxons(sp, taxon_map):
    sp = sp.replace("("," ")
    sp = sp.replace(")"," ")
    keywords = [kw.strip().lower() for kw in sp.split(" ") if kw.strip() != ""]
    
    strain_kws = set(['str.', 'substr.', 'strain', 'isolate'])
    tx_scores = []
    
    for k in taxon_map:
        for tx in taxon_map[k]:
            strains = []
            
            if tx.strain != None:        
                strain = tx.strain
                if strain.find("strain ") == 0:
                    strain = strain[len("strain "):]
                if strain.find("isolate ") == 0:
                    strain = strain[len("isolate "):]
                
                for st in strain.split("/"):
                    strains.extend([s.strip() for s in st.strip().split()])
            
            search_set = [s.lower() for s in tx.name.split(" ") + strains]
            
            score = 0
            for k in keywords:
                if k in search_set and k not in strain_kws:
                    score+=search_set.count(k)
                if k not in strain_kws:
                    score-=0.5
                if k in strain_kws and tx.strain != None:
                    score+=1
                    
            tx_scores.append((score, tx))
        
    return sorted(tx_scores, key=lambda (v, _): -v)
    

if __name__=='__main__':
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'development.ini'))
            
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        
        
        infile = open('data/speclist.txt', 'rb')
        
        node_id=None
        kingdom=None
        name=None        
        
        taxon_map = {}
        
        started=False
        for line in infile:
            if started:
                m = re.match(r"^[A-Z0-9]+\s+([A-Z])\s+([0-9]+):\s+N=([A-Za-z0-9\- ]*)(?:\((.*)\))?$", line)
                if m != None:
                    taxon = taxonomies.Taxonomy()
                    taxon.node_id = int(m.group(2))
                    taxon.kingdom = m.group(1)
                    taxon.name = m.group(3)
                    taxon.strain = m.group(4)
                    
                    t_set = taxon_map.get(taxon.name.lower(), set())
                    t_set.add(taxon)
                    taxon_map[taxon.name.lower()] = t_set
                    
                    DBSession.add(taxon)
            
            if line.find('_____') == 0:
                started=True
        infile.close()
        
        species = DBSession.query(taxonomies.Species).all()
        
        unmatched_species = []
        for sp in species:
            sp_name = sp.name.lower()
            if sp_name in taxon_map:
                tx_set = taxon_map[sp_name]
                if len(tx_set) == 1:
                    tx = list(tx_set)[0]
                    sp.taxon_id = tx.node_id
                    DBSession.add(sp)
                else:
                    unmatched_species.append(sp)
            else:
                unmatched_species.append(sp)
        
        print "%d unmatched species found " % (len(unmatched_species))
        
        ofile = open('output/unmatched_ranks', 'w')
        
        fixed_species = {}
        for sp in unmatched_species:
            ofile.write("%-90s" % (sp.name)) 
            similar = get_similar_taxons(sp.name, taxon_map)
            for i in xrange(0, 10):
                v, tx = similar[i]
                ofile.write("<%f, %d, %s, %s> -- " % (v, tx.node_id, tx.name, tx.strain))
            
            v, tx = similar[0]
            if v >= 1.0:
                fixed_species[sp.name] = (v, tx)
                sp.taxon_id = tx.node_id
                DBSession.add(sp)
            
            ofile.write("\n")
        ofile.close()
        
        ofile = open('output/inferred_assignments', 'w')
        for sp_name in fixed_species:
            v, tx = fixed_species[sp_name]
            if tx.strain != None:
                ofile.write("%90s -- %.4f %8d %s (strain %s)\n" % (sp_name, v, tx.node_id, tx.name, tx.strain))
            else:
                ofile.write("%90s -- %.4f %8d %s\n" % (sp_name, v, tx.node_id, tx.name))
        ofile.close()
        
        print "Fixed %d species assignments, remaining broken: %d" % (len(fixed_species), len(unmatched_species) - len(fixed_species))
        
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()