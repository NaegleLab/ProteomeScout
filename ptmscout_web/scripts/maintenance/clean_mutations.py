from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein, mutations
from paste.deploy.loadwsgi import appconfig
from ptmworker.helpers import upload_helpers
from ptmscout.config import settings
from geeneus import Proteome
import traceback
import sys, os

FLUSH_FREQ = 1000
if __name__ == "__main__":
    try:
        dbconfig = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'test.ini'))
        
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        prot_cnt = DBSession.query(protein.Protein).count()

        i = 0
        print "Getting protein data..."
        for prot in DBSession.query(protein.Protein):
            unique_mutations = {}
            mutation_dict = {}
            for m in prot.mutations:
                key = (m.mutationType, m.location, m.mutant, m.annotation)
                unique_mutations[key] = m.id
                mutation_dict[m.id] = m

            if len(unique_mutations) < len(prot.mutations):
                print "Found protein %d with revised mutations..." % ( prot.id )

#                print "Old:"
#                for m in prot.mutations:
#                    print m.location, m.original, m.mutant, m.annotation

                m_ids = set(unique_mutations.values())
                prot.mutations = []

                for m_id in mutation_dict:
                    if m_id in m_ids:
                        prot.mutations.append(mutation_dict[m_id])
                    else:
                        mutation_dict[m_id].protein_id = None
                        DBSession.delete(mutation_dict[m_id])

                DBSession.add(prot)

#                print "New:"
#                for m in prot.mutations:
#                    print m.location, m.original, m.mutant, m.annotation

            i+=1
            if i % FLUSH_FREQ == 0:
                DBSession.flush()
                print "%d / %d" % (i, prot_cnt)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
