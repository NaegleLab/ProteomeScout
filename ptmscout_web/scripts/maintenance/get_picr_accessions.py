from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein
from ptmworker.helpers import upload_helpers, picr_tools
import sys, os
import traceback

def get_picr_results(accessions, taxon_id):
    result = []
    for acc in accessions:
        result.extend( picr_tools.get_picr(acc, taxon_id) )
    return result

def get_primary_accessions(accessions):
    uniprot = []
    refseq = []
    for acc in accessions:
        if acc.type == 'uniprot':
            uniprot.append(acc.value)
        elif acc.type == 'refseq':
            refseq.append(acc.value)
    if uniprot != []:
        return uniprot

    return refseq

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

        for prot in proteins:
            paccessions = get_primary_accessions(prot.accessions)
            picr_accessions = get_picr_results(paccessions, prot.species.taxon_id)
            
            upload_helpers.create_accession_for_protein(prot, picr_accessions)
            
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
