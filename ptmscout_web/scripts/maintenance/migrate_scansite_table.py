from scripts.DB_init import DatabaseInitialization
from scripts import progressbar
from ptmscout.database import DBSession, modifications, protein
import os, sys
import traceback

YIELD_PER = 10000

def page_query(query_generator, start=0):
    offset = start
    while True:
        r = False
        for elem in query_generator(YIELD_PER, offset):
            r = True
            yield elem
            offset += YIELD_PER

        if not r: break

def page_peptides():
    def query_generator(limit, offset):
        sq = DBSession.query(modifications.Peptide.id).order_by(modifications.Peptide.id).limit(limit).offset(offset).subquery()
        return DBSession.query(modifications.Peptide).join(sq, modifications.Peptide.id == sq.c.id)
    return page_query( query_generator )

if __name__ == "__main__":
    try:
        settings = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', sys.argv[1])
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        total_peptides = DBSession.query(modifications.Peptide).count()
        i = 0

        pb = progressbar.ProgressBar(max_value = total_peptides)
        pb.start()


        j = 0
        for pep in page_peptides():
            for pred in pep.predictions:
                ps = protein.ProteinScansite()
                
                ps.percentile = pred.percentile
                ps.value = pred.value
                ps.source = pred.source
                ps.score = pred.score
                ps.site_pos = pep.site_pos
                
                if not pep.protein.hasPrediction(ps.source, ps.value, ps.site_pos):
                    pep.protein.scansite.append(ps)
                    j += 1
                
                pep.protein.saveProtein()
        
            i+=1
            pb.update(i)
            if i % 1000 == 0:
                DBSession.flush()
            
        pb.finish()

        sys.stderr.write("Wrote %d scansite annotations\n" % (j))
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
