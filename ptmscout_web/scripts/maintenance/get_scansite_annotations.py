from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, modifications, taxonomies
from paste.deploy.loadwsgi import appconfig
import time, datetime
from ptmworker import upload_helpers
import sys, os
from ptmscout.utils.decorators import rate_limit
import traceback

MOTIF_MAP = {
    'bos taurus':"MAMMALIAN",
    'bubalus bubalis':"MAMMALIAN",
    'callithrix jacchus':"MAMMALIAN",
    'canis familiaris':"MAMMALIAN",
    'cavia porcellus':"MAMMALIAN",
    'capra hircus':"MAMMALIAN",
    'chlorocebus aethiops':"MAMMALIAN",
    'cricetulus longicaudatus':"MAMMALIAN",
    'cricetulus griseus':"MAMMALIAN",
    'cricetus cricetus':"MAMMALIAN",
    'equus caballus':"MAMMALIAN",
    'homo sapiens':"MAMMALIAN",
    'mesocricetus auratus':"MAMMALIAN",
    'mus musculus':"MAMMALIAN",
    'mus musculus molossinus':"MAMMALIAN",
    'mustela putorius furo':"MAMMALIAN",
    'oryctolagus cuniculus':"MAMMALIAN",
    'ovis aries':"MAMMALIAN",
    'rattus norvegicus':"MAMMALIAN",
    'sus scrofa':"MAMMALIAN",
    'felis catus':"MAMMALIAN",
#    'Asterina pectinifera'
#    'Coturnix coturnix japonica'
#    'drosophila melanogaster'
#    'gallus gallus'
#    'Meleagris gallopavo'
#    'torpedo californica'
#    'xenopus laevis'
}

FLUSH_EVERY=100

@rate_limit(rate=10)
def insert_scansite_predictions(pep):
    species_name = pep.protein.species.name
    motif_class = MOTIF_MAP[species_name.lower()]

    predictions = upload_helpers.query_peptide_predictions(pep.pep_aligned, motif_class)

    if len(predictions) > 0:
        for p in predictions:
            pep.predictions.append(p)
        pep.scansite_date = datetime.datetime.now()
        DBSession.add(pep)



def needs_scansite(pep):
    species_name = pep.protein.species.name
    return species_name.lower() in MOTIF_MAP and pep.scansite_date==None


if __name__ == "__main__":
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'production.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        i = 0
        print "Getting peptide data..."

        peptides = DBSession.query(modifications.Peptide).all()
        for pep in peptides:
            if needs_scansite(pep):
                insert_scansite_predictions(pep)

            if i % FLUSH_EVERY == 0:
                print i, "/", len(peptides)
                DBSession.flush()
            i+=1
        
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
