from scripts.DB_init import DatabaseInitialization
from ptmscout.database import protein, modifications
import sys, os
import csv

def format_protein_accessions(accessions):
    return '; '.join([ acc.value for acc in sorted(accessions, key=lambda acc: acc.value) ])

def check_modtype_filter(mod, modtype_filter):
    if modtype_filter == None:
        return True

    modtype_filter = modtype_filter.lower()

    p = mod
    while(p is not None):
        if p.name.lower() == modtype_filter:
            return True
        for k in p.keywords:
            if k.keyword.lower() == modtype_filter:
                return True
        p = p.parent

    return False


def format_modifications(mods, modtype_filter):
    modlist = [ (modpep.peptide, modpep.modification) for ms in mods for modpep in ms.peptides ]
    modlist = [ (pep.site_pos, "%s-%s" % (pep.getName(), mod.name)) for pep, mod in modlist if check_modtype_filter(mod, modtype_filter) ]
    modlist = [ modstr for site_pos, modstr in sorted( list( set(modlist) ) ) ]
    n = len(modlist)
    return n, '; '.join(modlist)

def format_regions(regions):
    return '; '.join( [ "%s:%d-%d" % (r.label, r.start, r.stop) for r in sorted( regions, key=lambda r: r.start ) ] )

def format_mutations(mutations):
    return '; '.join( [ str(m) for m in sorted(mutations, key=lambda m: m.location) ] )

def format_GO_terms(prot):
    return '; '.join( [ "%s-%s" % (goe.GO_term.GO, goe.GO_term.term) for goe in sorted(prot.GO_terms, key=lambda term: term.GO_term.GO) ] )

def format_scansite(mods):
    plist = [ (modpep.peptide.site_pos, pred.source, pred.value, pred.percentile) for ms in mods for modpep in ms.peptides for pred in modpep.peptide.predictions ]
    return '; '.join( [ '%d %s %s (%.2f)' % (pos, source, value, percentile) for pos, source, value, percentile in sorted( list( set( plist ) ) ) ])

def check_species_filter(f, p):
    if f == None:
        return True
    f = f.lower()

    if f == p.species.name.lower():
        return True

    t = p.species.taxon
    while t != None:
        if f == t.name.lower():
            return True

        t = t.parent

    return False

if __name__=='__main__':
    species_filter = None
    modtype_filter = None

    for i in xrange(1, len(sys.argv)):
        if sys.argv[i] == '--species':
            species_filter = sys.argv[i+1]
        if sys.argv[i] == '--modification':
            modtype_filter = sys.argv[i+1]

    config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
    DatabaseInitialization.setUpClass(config_options)
    db = DatabaseInitialization()
    db.setUp()

    prot_cnt = db.session.query(protein.Protein).count()
    all_proteins = db.session.query(protein.Protein)

    cw = csv.writer(sys.stdout, dialect='excel-tab')

    cw.writerow(['accessions', 'acc_gene', 'locus', 'protein_name', 
                    'species', 'sequence', 'modifications', 'domains', 
                    'mutations', 'scansite_predictions', 'GO_terms'])

    i = 0
    j = 0
    k = 0
    for p in all_proteins:
        if check_species_filter(species_filter, p):
            mods = modifications.getMeasuredPeptidesByProtein(p.id)
            n, fmods = format_modifications(mods, modtype_filter)
            if n > 0:
                k+=n
                row = []
                row.append( format_protein_accessions(p.accessions) )
                row.append( p.acc_gene )
                row.append( p.locus )
                row.append( p.name )
                row.append( p.species.name )
                row.append( p.sequence )
                row.append( fmods )
                row.append( format_regions(p.domains) )
                row.append( format_mutations(p.mutations) )
                row.append( format_scansite(mods) )
                row.append( format_GO_terms(p) )
                cw.writerow(row)
                j+=1

        i+=1
        if i % 100 == 0:
            sys.stderr.write( "%d / %d\n" % (i, prot_cnt) )
    
    sys.stderr.write( 'Exported %d proteins, with %d unique modifications' % (j, k))
    db.close()
