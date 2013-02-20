from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein, gene_expression
import sys, os
import traceback
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

        protein_map = {}

        queries = []
        for prot in proteins:
            paccessions = [ acc.value.lower() for acc in prot.accessions ]
            probesets = gene_expression.getExpressionProbeSetsForProtein(paccessions, prot.species_id)
            
            added = 0
            for probeset in probesets:
                if probeset not in prot.expression_probes:
                    prot.expression_probes.append(probeset)
                    added+=1
            
            if added > 0:
                print "Added %d / %d probesets" % (added, len(probesets))
                DBSession.add(prot)

            i+=1
            if i % FLUSH_EVERY == 0:
                DBSession.flush()
                print i, "/", len(proteins)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
