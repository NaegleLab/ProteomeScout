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
PICKLEFILE = './mutations/mutationsUpdate.pickle'
FILEBASE = './mutations/accList_'
if __name__ == "__main__":
    try:
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        manager = Proteome.ProteinManager(settings.adminEmail, cache=True, retry=2) # but remember to purge!

	#prot_map = pickle.load(open(PICKLEFILE, "r"))
	#	pcnt = len(prot_map)
        
        i = 0
	fileNums = [1]
	queryBad = []
        for batch in fileNums:
            i+=1
	    # read in all of the accssions 
	    fileName = FILEBASE+str(batch)
	    with open (fileName, 'r') as f:
		batch_query = [line.rstrip() for line in f]
	    print "there are %d accessions in %s"%(len(batch_query), fileName)
	    try:
		    manager.batch_get_protein_name(batch_query)
	    except Exception, e2: 
		    print "Error in get_protein_name with error %s, attempting recovery"%(e2)
		    for query in batch_query:
			try:
				seq = manager.get_protein_sequence(query)
				if seq is None:
					print "Error getting sequence %s with error: %s, removing from protein map"%(query, e2)
					queryBad.append(query)
					protein_map.pop(query, None)
			except Exception, e3:
				print "Error getting sequence %s with error: %s, removing from protein map"%(query, e3)
				queryBad.append(query)
				batch_query.remove(query)
				protein_map.pop(query, None)
			else:
				print "Success for %s"%(query)
	    else: 
		print "No errors on that try catch for batch %s" %(fileName)
	
	print "Total bad queries:\n"
	print queryBad
        print "Building variant table..."
        i = 0
        for acc in batch_query:
	    prots = protein.getProteinsByAccession(acc)
            #prot = prot_map[acc]
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
                            print "Loaded mutation does not match protein sequence %s (%d %s) %s -> %s" % (acc, new_mutation.location, prot_seq[new_mutation.location-1], new_mutation.original, new_mutation.mutant)
                        elif not prot.hasMutation(new_mutation):
                            j+=1
                            DBSession.add(new_mutation)

                print "Added %d mutations for protein %d %s"  % (j, prot.id, acc)
            else:
                print "Warning: protein sequence mismatch for %d with %s" % (prot.id, acc)

            i+=1
            if i % FLUSH_FREQ == 0:
                DBSession.flush()

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
	print "Completed access list %s"%(fileName)
        dbinit.tearDown()
