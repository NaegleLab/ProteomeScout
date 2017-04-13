from scripts.DB_init import DatabaseInitialization
from ptmscout.database import protein, modifications, experiment
import sys, os
import csv
from scripts.progressbar import ProgressBar
from ptmscout.utils.export_proteins import *
from ptmscout.utils.downloadutils import experiment_metadata_to_tsv, zip_package

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

    output_fn_root = sys.argv[-1]

    config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
    DatabaseInitialization.setUpClass(config_options)
    db = DatabaseInitialization()
    db.setUp()

    prot_cnt = db.session.query(protein.Protein).count()
    all_proteins = db.session.query(protein.Protein.id)

    if not os.path.exists(output_fn_root):
        os.makedirs(output_fn_root)
    data_filename = os.path.join(output_fn_root, "data.tsv")
    metadata_filename = os.path.join(output_fn_root, "citations.tsv")
    zip_filename = output_fn_root + ".zip"

    experiment_list = set()

    with open(data_filename, 'w') as dfile:
        cw = csv.writer(dfile, dialect='excel-tab')

        cw.writerow(['protein_id', 'accessions', 'acc_gene', 'locus', 'protein_name',\
                        'species', 'sequence', 'modifications', 'evidence', \
                        'pfam_domains', 'uniprot_domains',\
                        'kinase_loops', 'macro_molecular',\
                        'topological', 'structure',\
                        'scansite_predictions', 'GO_terms',\
                        'mutations','mutation_annotations'])

        pb = ProgressBar(max_value = prot_cnt)
        pb.start()
        i = 0
        j = 0
        k = 0
        for p_id in all_proteins:
	    p = db.session.query(protein.Protein).get(p_id)
            if check_species_filter(species_filter, p):
                mods = modifications.getMeasuredPeptidesByProtein(p.id)
                qaccs = get_query_accessions(mods)
                n, fmods, fexps, fexp_id_set = format_modifications(mods, modtype_filter)
                experiment_list |= fexp_id_set
                if n > 0:
                    k+=n
                    row = []
                    row.append( p.id )
                    row.append( format_protein_accessions(p.accessions, qaccs) )
                    row.append( p.acc_gene )
                    row.append( p.locus )
                    row.append( p.name )
                    row.append( p.species.name )
                    row.append( p.sequence )
                    row.append( fmods )
                    row.append( fexps )

                    uniprot_domains = filter_regions(p.regions, set(['domain']))
                    kinase_loops = filter_regions(p.regions, set(['Activation Loop']))
                    macromolecular = filter_regions(p.regions, set([ 'zinc finger region', 'intramembrane region', 'coiled-coil region', 'transmembrane region' ]))
                    topological = filter_regions(p.regions, set(['topological domain']))
                    structure = filter_regions(p.regions, set(['helix', 'turn', 'strand']))

                    row.append( format_domains(p.domains) )
                    row.append( format_domains(uniprot_domains) )
                    row.append( format_domains(kinase_loops) )
                    row.append( format_regions(macromolecular) )
                    row.append( format_domains(topological) )
                    row.append( format_regions(structure) )

                    row.append( format_scansite(mods) )
                    row.append( format_GO_terms(p) )
                    row.append( format_mutations(p.mutations) )
                    row.append( format_mutation_annotations(p.mutations) )
                    cw.writerow(row)
                    j+=1

            i+=1
            if i % 10 == 0:
                pb.update(i)

    pb.finish()
    sys.stderr.write( 'Exporting experiment metadata...' )

    experiments = [ experiment.getExperimentById(exp_id) for exp_id in experiment_list ]
    experiment_metadata_to_tsv(experiments, metadata_filename)

    zip_package([data_filename, metadata_filename], zip_filename)

    sys.stderr.write( 'Exported %d proteins, with %d unique modifications' % (j, k))
    db.close()
