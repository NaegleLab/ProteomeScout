from DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications, taxonomies
from paste.deploy.loadwsgi import appconfig
import time
from ptmworker import upload_helpers
import sys
from ptmscout.utils.decorators import rate_limit

MOTIF_MAP = {
    'bos taurus':"MAMMALIAN",
    'capra hircus':"MAMMALIAN",
    'cricetulus griseus':"MAMMALIAN",
    'homo sapiens':"MAMMALIAN",
    'mus musculus':"MAMMALIAN",
    'oryctolagus cuniculus':"MAMMALIAN",
    'ovis aries':"MAMMALIAN",
    'rattus norvegicus':"MAMMALIAN",
    'sus scrofa':"MAMMALIAN",
    'Felis catus':"MAMMALIAN"
#    'drosophila melanogaster'
#    'gallus gallus'
#    'xenopus laevis'
}

@rate_limit(rate=10)
def insert_scansite_predictions(pep):
    species_name = pep.protein.species.name
    motif_class = MOTIF_MAP[species_name]

    predictions = upload_helpers.query_peptide_predictions(pep.pep_aligned, motif_class)

    if len(predictions) > 0:
        pep.predictions.extend(predictions)
        DBSession.add(pep)



def needs_scansite(pep):
    species_name = pep.protein.species.name
    return species_name in MOTIF_MAP and len(pep.predictions) == 0


if __name__ == "__main__":
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'production.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        i = 0
        for pep in DBSession.query(modifications.Peptide):
            if needs_scansite(pep):
                insert_scansite_predictions(pep)

            i+=1
            if i % FLUSH_EVERY == 0:
                print i
                DBSession.flush()
        
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
