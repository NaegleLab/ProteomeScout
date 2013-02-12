from Bio import SeqIO
import sys, os
from scripts.DB_init import DatabaseInitialization
from ptmscout.database import DBSession, uniprot
import re
import traceback

# sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) GN=FV3-001R PE=4 SV=1
# sp|P48347-2|14310_ARATH Isoform 2 of 14-3-3-like protein GF14 epsilon OS=Arabidopsis thaliana GN=GRF10

def get_last_param(param_string):
    # OS, SV, PE, GN
    items = [0,0,0,0]
    items[0] = param_string.rfind('OS=')
    items[1] = param_string.rfind('SV=')
    items[2] = param_string.rfind('PE=')
    items[3] = param_string.rfind('GN=')

    last = sorted(items)[-1]

    if last == -1:
        return None

    p = param_string[last:last+2]
    v = param_string[last+3:]
    rest = param_string[:last]
    return p, v.strip(), rest


def parse_params(name):
    params = {}

    try:
        while(1):
            p, v, remaining = get_last_param(name)
            name = remaining
            params[p] = v
    except TypeError:
        return name.strip(), params



def parse_record_description(desc):
    m = re.match(r"^sp\|([A-Za-z0-9\_\-]+)\|([A-Za-z0-9\_\-]+) (.*)$", desc)

    if m:
        accession = m.group(1)
        locus = m.group(2)

        name, params = parse_params(m.group(3))
        species = params['OS']

        return name, accession, locus, species


def load_swissprot(filename):
    handle = open(filename,'rU')

    print "Processing '%s'..." % (filename)
    i = 0
    for record in SeqIO.parse(handle, "fasta"):
        result = parse_record_description(record.description)
        
        if result:
            name, accession, locus, species = result
            sp_record = uniprot.SwissprotRecord(name, accession, locus, species, record.seq)
            DBSession.add(sp_record)
            i+=1
        else:
            print "Failed to create record for entry: %s " % (record.description)


    print "Created %d records" % ( i )
    handle.close()

if __name__ == "__main__":
    swissprot_fn = sys.argv[1]
    isoform_fn = sys.argv[2]
    database = sys.argv[3]

    if database == 'test':
        dbconfig = os.path.join('data', 'ptmscout', 'ptmscout_web', 'test.ini')
    elif database == 'production':
        dbconfig = os.path.join('data', 'ptmscout', 'ptmscout_web', 'production.ini')


    try:
        DatabaseInitialization.setUpClass(dbconfig)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        load_swissprot(swissprot_fn)
        load_swissprot(isoform_fn)

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
