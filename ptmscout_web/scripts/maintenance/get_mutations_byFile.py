from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, Protein, taxonomies, protein, mutations
from ptmscout.config import settings
from geeneus import Proteome
import traceback
import os
import pickle


DEBUG=False
REPORT_EVERY = 250
BATCH_SIZE = 500
FLUSH_FREQ = 1000
LOGFILE = './mutations/logFile'
FILEBASE = './mutations/accList_'

def get_mutations_for_accList(accFile, logFile):
    flag = 0
    fLog = open(logFile, 'a')
    try:

	manager = Proteome.ProteinManager(settings.adminEmail, cache=True, retry=2) # but remember to purge!
        i = 0
	queryBad = []
	i+=1
	with open (accFile, 'r') as f:
		batch_query = [line.rstrip() for line in f]
		fLog.write("there are %d accessions in %s\n"%(len(batch_query), accFile))
	try:
		manager.batch_get_protein_name(batch_query)
	except Exception, e2: 
		fLog.write("Error in get_protein_name with error %s, attempting recovery\n"%(e2))
		for query in batch_query:
			try:
				seq = manager.get_protein_sequence(query)
				if seq is None:
					fLog.write("Error getting sequence %s with error: %s, removing from protein map\n"%(query, e2))
					queryBad.append(query)
					batch_query.remove(query)
			except Exception, e3:
				fLog.write("Error getting sequence %s with error: %s, removing from protein map\n"%(query, e3))
				queryBad.append(query)
				batch_query.remove(query)
			else:
				fLog.write("Success for %s\n"%(query))
	else: 
		fLog.write("No errors on that try catch for batch %s\n" %(accFile))
	
	fLog.write("Total bad queries: for %s is %d \n"%(accFile, len(queryBad)))
	for q in queryBad:
		fLog.write("%s, "%(q))
	fLog.write("\n")
        i = 0
        for acc in batch_query:
	    prots = protein.getProteinsByAccession(acc)
            prot = prots[0]
            prot_seq = manager.get_protein_sequence(acc)
            if prot_seq and prot.sequence == prot_seq.upper():
                prot_seq = prot_seq.upper()
                variantTuple = manager.get_variants(acc)
                
                j = 0
                for mutantDict in variantTuple:
                    # for now we're only working with single mutants, but could expand
                    # to double mutants in the future...
                    if(mutantDict['type'] == "Substitution (single)"):
                        new_mutation = mutations.Mutation(mutantDict['type'],
                                mutantDict['location'], mutantDict['original'],
                                mutantDict['mutant'], acc,
                                mutantDict['notes'], prot.id)
                        if not new_mutation.consistent(prot_seq):
                            fLog.write("Loaded mutation does not match protein sequence %s (%d %s) %s -> %s\n" % (acc, new_mutation.location, prot_seq[new_mutation.location-1], new_mutation.original, new_mutation.mutant))
                        elif not prot.hasMutation(new_mutation):
                            j+=1
                            DBSession.add(new_mutation)

                fLog.write("Added %d mutations for protein %d %s\n"  % (j, prot.id, acc))
            else:
                fLog.write("Warning: protein sequence mismatch for %d with %s\n" % (prot.id, acc))

            i+=1
            if i % FLUSH_FREQ == 0:
                DBSession.flush()

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        fLog.write("Rolling back database changes... for %s\n"%(accFile))
        dbinit.rollback()
	flag = -1
    else:
        print "Finalizing DB changes"
	flag = 1
	print "Completed access list %s\n"%(accFile)
        dbinit.tearDown()
    fLog.close()
    return flag, queryBad

if __name__ == "__main__":
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
	
	j = 0
	accFile = FILEBASE+str(j)
	flag, queryBad = get_mutations_for_accList(accFile, LOGFILE)
	if(flag==1): 
		print "Success for %s"%(accFile)
