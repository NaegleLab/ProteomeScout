from scripts import progressbar
from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, protein, mutations
from ptmscout.config import settings
import traceback
import os, sys
from ptmworker.protein_tasks import load_new_protein, get_proteins_by_accession

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

if __name__ == "__main__":
    try:
        dbconfig = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        accessions = []
        for prot in DBSession.query(protein.Protein):
            pacc = get_primary_accessions(prot)
            if len(pacc) > 0:
                acc = pacc[0]
                accessions.append(acc)
    
        sys.stderr.write( "Querying for %d accessions\n" % (len(accessions)) )
        pb = progressbar.ProgressBar()

        def notify_callback(i, total, errors):
            pb.update(i)

        def start_callback(total):
            pb.max_value = total
            pb.start()

        protein_map = get_proteins_by_accession(dict((a,None) for a in accessions), start_callback, notify_callback)
        pb.finish()

        sys.stderr.write( "Loading queried proteins\n" )
        pb = progressbar.ProgressBar(max_value = len(protein_map))
        pb.start()

        i = 0
        mismatches = 0
        loaded = 0
        errors = {}
        for acc in protein_map:
            protein_record = protein_map[acc]
            prot = protein.getProteinBySequence(protein_record.sequence, protein_record.species)
            if prot == None:
                mismatches+=1
            else:
                try:
                    load_new_protein(acc, protein_record)
                    loaded+=1
                except Exception, e:
                    errors[acc] = e

            i+=1
            pb.update(i)
        pb.finish()

        sys.stderr.write("Reloaded info for %d proteins (%d mismatches, %d errors)" % (loaded, mismatches, len(errors)))

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
