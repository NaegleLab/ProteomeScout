from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, Protein, taxonomies, protein, mutations
from ptmscout.config import settings
from geeneus import Proteome
import traceback
import os

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
    try:
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        manager = Proteome.ProteinManager(settings.adminEmail, cache=True, retry=2) # but remember to purge!
        pcnt = DBSession.query(Protein).join(Protein.species).filter(taxonomies.Species.name=='homo sapiens').count()
        q = DBSession.query(Protein).join(Protein.species).filter(taxonomies.Species.name=='homo sapiens').count()

        
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
            print "querying for accession batch... %d / %d" % (i, len(batch_queries))
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
				protein_map.pop(query, None)
			else:
				print "Success for %s"%(query)
	    else: 
		print "No errors on that try catch for batch %d / %d" %(i, len(batch_queries))
	
		print "Total bad queries:\n"
		print queryBad
        print "Building variant table..."

        print "Building variant table..."
        i = 0
        for acc in protein_map:
            prot = protein_map[acc]
            
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
        dbinit.tearDown()
