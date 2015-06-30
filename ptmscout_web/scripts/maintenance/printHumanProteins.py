from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, Protein, taxonomies, protein, mutations
from ptmscout.config import settings
from geeneus import Proteome
import traceback
import os
import pickle

def get_primary_accessions(prot):
    uniprot_accessions = []
    refseq_accessions = []
    genbank_accessions = []
    ipi_accessions = []
    for acc in prot.accessions:
        tp = acc.type.lower()
        if tp == 'uniprot' or tp == 'swissprot':
            uniprot_accessions.append(acc.value)
        if tp == 'refseq':
            refseq_accessions.append(acc.value)
        if tp == 'genbank':
            genbank_accessions.append(acc.value)
        if tp == 'ipi':
            genbank_accessions.append(acc.value)

    primary_accessions = sorted(uniprot_accessions, key=lambda item: len(item)) + refseq_accessions
    secondary_accessions = genbank_accessions + ipi_accessions

    return primary_accessions if len(primary_accessions) > 0 else secondary_accessions

DEBUG=False
REPORT_EVERY = 250
BATCH_SIZE = 500
FLUSH_FREQ = 1000
if __name__ == "__main__":
# setup the pickle file for the session and print accessions to files, so that individual files can be run for mutations, this can be done in parallel and then we'll know which protein sets have been updated.
    try:
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        manager = Proteome.ProteinManager(settings.adminEmail, cache=True, retry=2) # but remember to purge!
	pcnt = DBSession.query(Protein).join(Protein.species).filter(taxonomies.Species.name=='homo sapiens').count()
	q = DBSession.query(Protein).join(Protein.species).filter(taxonomies.Species.name=='homo sapiens').count()
	pickleFile = './mutations/mutationsUpdate.pickle'
	fileBase = './mutations/accList'

	
        
        protein_map = {}
        batch_queries = [[]]

        i = 0
        print "Getting protein data..."
        for prot in DBSession.query(Protein).join(Protein.species).filter(taxonomies.Species.name=='homo sapiens'):
            accessions = get_primary_accessions(prot)
            if len(accessions) == 0:
                print "Warning: No prefered or secondary accessions for record %d [ %s ]" % (prot.id, ', '.join([acc.value for acc in prot.accessions]))
            elif DEBUG:
                print "Adding %d accessions for protein %d %s" % ( len(accessions), prot.id, ', '.join([acc for acc in accessions]) )

            for acc in accessions:
                protein_map[acc] = prot
                batch_queries[-1].append(acc)

                if len(batch_queries[-1]) == BATCH_SIZE:
                    batch_queries.append([])

            i+=1
            if i % REPORT_EVERY == 0:
                print "Protein %d / %d" % (i, pcnt)

        i = 0
        for batch_query in batch_queries:
		i+=1
		fileOut = fileBase+'_'+str(i)
		print "Writing file %s"%(fileOut)
		f = open(fileOut, 'w')
		for acc in batch_query:
			f.write(acc+"\n")
		f.close()	
	    
	pickle.dump(protein_map, open(pickleFile, "w"))	
    except Exception, e:
	print "error %s"%(e)
