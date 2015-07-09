from scripts import progressbar
from ptmworker.helpers import dbsnp_tools
from ptmscout.database import mutations
from scripts.DB_init import DatabaseInitialization
from ptmscout.config import settings
from ptmscout.database import mutations
import traceback
import os, sys

QUERY_SIZE = 100

def query_mutation_info(snps, session):
    result = dbsnp_tools.get_variant_classification(snps.keys())
    for rsid in snps:
        if rsid in result:
            snps[rsid].clinical = result[rsid]
            session.add(snps[rsid])
    return len(result)

if __name__ == "__main__":
    try:
	
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        cnt = dbinit.session.query(mutations.Mutation).count()
        pb = progressbar.ProgressBar(max_value=cnt)
        pb.start()

        dbsnps = {}

        total = 0
        i = 0
        for m in dbinit.session.query(mutations.Mutation):
            dbsnp_id = m.getDBSNPId()

            if dbsnp_id != None:
                dbsnps[dbsnp_id] = m

            if len(dbsnps) == QUERY_SIZE:
                total += query_mutation_info(dbsnps, dbinit.session)
                dbsnps = {}
                dbinit.session.flush()

            i+=1
            if i % 100 == 0:
                pb.update(i)

        if len(dbsnps) > 0:
            total += query_mutation_info(dbsnps, dbinit.session)

        sys.stderr.write("Queried: %d dbSNP entries\n" % (total))
        pb.finish()

        dbinit.session.flush()
    except Exception, e:
        traceback.print_exc()
        dbinit.rollback()
    else:
        dbinit.tearDown()
