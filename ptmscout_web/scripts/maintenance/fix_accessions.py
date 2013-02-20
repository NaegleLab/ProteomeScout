from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
import sys, os
import traceback
import re

FLUSH_EVERY = 100

if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        i = 0
        print "Getting protein data..."

        if len(sys.argv) > 1:
            query_ids = [int(pid) for pid in sys.argv[1:]]
            proteins = DBSession.query(protein.Protein).filter(protein.Protein.id.in_(query_ids)).all()
        else:
            proteins = DBSession.query(protein.Protein).all()

        type_map = {'ref':'refseq', 'gb':'genbank', 'emb': 'EMBL'}

        print "Fixing accessions..."
        for prot in proteins:
            remove = []
            for acc in prot.accessions:
                m = re.match('^([a-zA-Z]+)\|(.*)\|$', acc.value)
                if m:
                    tp = m.group(1)
                    val = m.group(2)
                    if prot.hasAccession(val):
                        remove.append(acc)
                        continue

                    if tp in type_map:
                        print "Changed accession: ( %s, %s )  ==>  ( %s, %s )" % ( str(acc.type), str(acc.value), str(type_map[tp]), str(val) )
                        acc.type = type_map[tp]
                        acc.value = val
                    else:
                        remove.append(acc)

            for acc in remove:
                print "Removing duplicate accession ( %s, %s )" % ( str(acc.type), str(acc.value) )
                prot.accessions.remove(acc)

            prot.saveProtein()

            i+=1
            if i % FLUSH_EVERY == 0:
                DBSession.flush()

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
