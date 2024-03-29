from ptmworker.helpers import upload_helpers
from ptmscout.config import settings, strings
from ptmscout.database import experiment, modifications, user, modifications, protein
import celery
from ptmworker import notify_tasks, protein_tasks
from ptmscout.utils import export_proteins, downloadutils
import csv, os, random

NOTIFY_INTERVAL = 5

def annotate_experiment(user, exp, header, rows, job_id):
    notify_tasks.set_job_stage.apply_async((job_id, 'annotating', len(rows)))

    header += [ 'scansite_bind', 'scansite_kinase', 'nearby_modifications',\
            'nearby_mutations', 'nearby_mutation_annotations', \
                    'site_pfam_domains', 'site_uniprot_domains',\
                    'site_kinase_loop', 'site_macro_molecular',\
                    'site_topological', 'site_structure',\
                    'protein_pfam_domains', 'protein_uniprot_domains',\
                    'protein_GO_BP', 'protein_GO_CC', 'protein_GO_MF' ]
    protein_mods = {}
    
    for ms in exp.measurements:
        if ms.protein_id not in protein_mods:
            protein_mods[ms.protein_id] = modifications.getMeasuredPeptidesByProtein(ms.protein_id, user)
    
    ms_map = {}
    for ms in exp.measurements:
        ms_map[ms.id] = ms
        
    i = 0
    mx_val = len(rows)
    for row in rows:
        ms = ms_map[row[0]]
        prot = ms.protein
        
        min_range = ms.peptides[0].peptide.site_pos - 7
        max_range = ms.peptides[-1].peptide.site_pos + 7
        
        nearby_modifications = set()
        for ms2 in protein_mods[ms.protein_id]:
            for modpep in ms2.peptides:
                site_type = modpep.peptide.site_type
                site_pos = modpep.peptide.site_pos
                mod_name = modpep.modification.name
                
                if min_range <= site_pos and site_pos <= max_range:
                    nearby_modifications.add((site_pos, site_type, mod_name))
                    
        nearby_modifications = [ "%s%d: %s" % (site_type, site_pos, mod_name) for site_pos, site_type, mod_name in sorted(list(nearby_modifications)) ]
        nearby_mutations = [ mutation for mutation in sorted(prot.mutations, key=lambda item: item.location) if min_range < mutation.location and mutation.location < max_range ]
        
        sep = settings.mod_separator_character + ' '
        
        scansite_kinase = []
        scansite_bind = []
        for modpep in ms.peptides:
            for pp in modpep.peptide.predictions:
                if pp.source=='scansite_kinase':
                    scansite_kinase.append( "%s (%.2f)" % ( pp.value, pp.percentile ))
                if pp.source=='scansite_bind':
                    scansite_bind.append( "%s (%.2f)" % ( pp.value, pp.percentile ))

        pfam_sites = export_proteins.filter_sites(ms, prot.domains)
        domain_sites = export_proteins.filter_site_regions(ms, prot.regions, set(['domain']))
        kinase_sites = export_proteins.filter_site_regions(ms, prot.regions, set(['Activation Loop']))
        macromolecular_sites = export_proteins.filter_site_regions(ms, prot.regions, set([ 'zinc finger region', 'intramembrane region', 'coiled-coil region', 'transmembrane region' ]))
        topological_sites = export_proteins.filter_site_regions(ms, prot.regions, set(['topological domain']))
        site_structure = export_proteins.filter_site_regions(ms, prot.regions, set(['helix', 'turn', 'strand']))

        protein_uniprot_domains = export_proteins.filter_regions(prot.regions, set(['domain']))

        protein_GO = { 'P':set(), 'F':set(), 'C':set() }

        for ge in prot.GO_terms:
            term = ge.GO_term
            protein_GO[term.aspect].add(term.GO)

        row.append( sep.join(scansite_bind) )
        row.append( sep.join(scansite_kinase) )

        row.append( sep.join(nearby_modifications) )
        row.append( export_proteins.format_mutations( nearby_mutations ) )
        row.append( export_proteins.format_mutation_annotations( nearby_mutations ) )

        row.append( export_proteins.format_domains( pfam_sites ) )
        row.append( export_proteins.format_domains( domain_sites ) )
        row.append( export_proteins.format_domains( kinase_sites ) )
        row.append( export_proteins.format_regions( macromolecular_sites ) )
        row.append( export_proteins.format_domains( topological_sites ) )
        row.append( export_proteins.format_regions( site_structure ) )

        row.append( export_proteins.format_domains(prot.domains) )
        row.append( export_proteins.format_domains(protein_uniprot_domains) )
        
        row.append( sep.join(list(protein_GO['P'])) )
        row.append( sep.join(list(protein_GO['C'])) )
        row.append( sep.join(list(protein_GO['F'])) )

        i+=1
        if i % NOTIFY_INTERVAL == 0:
            notify_tasks.set_job_progress.apply_async((job_id, i, mx_val))

    return header, rows
        
def get_experiment_header(exp):
    header = ['MS_id', 'query_accession', 'gene', 'locus', 'protein_name', 'species', 'peptide', 'mod_sites', 'gene_site', 'aligned_peptides', 'modification_types']
    
    data_labels = set()
    for ms in exp.measurements:
        for d in ms.data:
            data_labels.add((d.run,d.type,d.units,d.label))
    
    def float_last_term(r,dt,u,l):
        try:
            l = float(l)
        except:
            pass
        
        return (r,dt,u,l)
    
    data_labels = [ "%s:%s:%s:%s" % (r,dt,u,str(l)) for r,dt,u,l in sorted(list(data_labels), key=lambda item: float_last_term(*item)) ]
    header += data_labels
    
    return header, data_labels

def get_experiment_data(exp, data_labels):
    rows = []
    for ms in exp.measurements:
        mod_sites = '; '.join([modpep.peptide.getName() for modpep in ms.peptides])
        aligned_peptides = '; '.join([modpep.peptide.pep_aligned for modpep in ms.peptides])
        modification_types = '; '.join([modpep.modification.name for modpep in ms.peptides])
        
        gene_sites = [ms.protein.getGeneName()] + [modpep.peptide.getName() for modpep in ms.peptides]
        row = [ms.id, ms.query_accession, ms.protein.acc_gene, ms.protein.locus, ms.protein.name, ms.protein.species.name, ms.peptide, mod_sites, '_'.join(gene_sites), aligned_peptides, modification_types]
        
        ms_data = {}
        for d in ms.data:
            formatted_label = "%s:%s:%s:%s" % (d.run, d.type, d.units, d.label)
            ms_data[formatted_label] = d.value
            
        for dl in data_labels:
            row.append( ms_data[dl] )
            
        rows.append(row)
            
    return rows


@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.dynamic_transaction_task
def run_experiment_export_job(annotate, export_id, exp_id, user_id, job_id):
    notify_tasks.set_job_status.apply_async((job_id, 'started'))
    notify_tasks.set_job_stage.apply_async((job_id, 'exporting', 0))

    exp_filename = 'experiment.%d.%d.%d.tsv' % (exp_id, user_id, export_id)
    exp_path = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, exp_filename)

    usr = user.getUserById(user_id)
    exp = experiment.getExperimentById(exp_id, usr)

    header, data_labels = get_experiment_header(exp)
    rows = get_experiment_data(exp, data_labels)

    if annotate:
        header, rows = annotate_experiment(usr, exp, header, rows, job_id)

    with open(exp_path, 'w') as export_file:
        cw = csv.writer(export_file, dialect='excel-tab')

        cw.writerow(header)
        for row in rows:
            cw.writerow(row)

    finalize_task = notify_tasks.finalize_experiment_export_job.s()
    return finalize_task, (job_id,), None

@celery.task
@upload_helpers.notify_job_failed
def annotate_proteins(protein_result, accessions, batch_id, exp_id, user_id, job_id):
    protein_map, protein_id_map = protein_result

    usr = user.getUserById(user_id)
    notify_tasks.set_job_stage.apply_async((job_id, 'annotate', len(protein_map)))
    data_filename = "batch.data.%s.%d.tsv" % (batch_id, user_id)
    metadata_filename = "batch.metadata.%s.%d.tsv" % (batch_id, user_id)
    zip_filename =  "batch.%s.%d.zip" % (batch_id, user_id)

    data_filepath = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, data_filename)
    metadata_filepath = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, metadata_filename)
    zip_filepath = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, zip_filename)

    header = ['protein_id', 'query_accession', 'other_accessions', 'acc_gene', 'locus', 'protein_name',\
                    'species', 'sequence', 'modifications', 'evidence',\
                    'pfam_domains', 'uniprot_domains',\
                    'kinase_loops', 'macro_molecular',\
                    'topological', 'structure',\
                    'mutations', 'mutation_annotations', 'scansite_predictions', 'GO_terms']
    rows = []
    success = 0
    errors = 0

    experiment_list = set()

    i = 0
    for acc in accessions:
        if acc in protein_map:
            pr = protein_map[acc]
            p = protein.getProteinBySequence( pr.sequence, pr.species )
            mods = modifications.getMeasuredPeptidesByProtein(p.id, usr)

            qaccs = export_proteins.get_query_accessions(mods)
            n, fmods, fexps, exp_list = export_proteins.format_modifications(mods, None)
            experiment_list |= exp_list

            row = []
            row.append( p.id )
            row.append( acc )
            row.append( export_proteins.format_protein_accessions(p.accessions, qaccs) )
            row.append( p.acc_gene )
            row.append( p.locus )
            row.append( p.name )
            row.append( p.species.name )
            row.append( p.sequence )
            row.append( fmods )
            row.append( fexps )

            uniprot_domains = export_proteins.filter_regions(p.regions, set([ 'domain' ]))
            kinase_loops = export_proteins.filter_regions(p.regions, set([ 'Activation Loop' ]))
            macromolecular = export_proteins.filter_regions(p.regions, set([ 'zinc finger region', 'intramembrane region', 'coiled-coil region', 'transmembrane region' ]))
            topological = export_proteins.filter_regions(p.regions, set([ 'topological domain' ]))
            structure = export_proteins.filter_regions(p.regions, set([ 'helix', 'turn', 'strand' ]))

            row.append( export_proteins.format_domains(p.domains) )
            row.append( export_proteins.format_domains(uniprot_domains) )
            row.append( export_proteins.format_domains(kinase_loops) )
            row.append( export_proteins.format_regions(macromolecular) )
            row.append( export_proteins.format_domains(topological) )
            row.append( export_proteins.format_regions(structure) )

            row.append( export_proteins.format_mutations(p.mutations) )
            row.append( export_proteins.format_mutation_annotations(p.mutations) )
            row.append( export_proteins.format_scansite(mods) )
            row.append( export_proteins.format_GO_terms(p) )
            rows.append( row )
            success += 1
        else:
            errors_for_acc = [ e.message for e in experiment.errorsForAccession(exp_id, acc) ]
            rows.append(['%d ERRORS: %s' % ( len(errors_for_acc), '; '.join(errors_for_acc) ), acc])
            errors+=1

        i+=1
        if i % NOTIFY_INTERVAL == 0:
            notify_tasks.set_job_progress.apply_async((job_id, i, len(protein_map)))

    with open(data_filepath, 'w') as bfile:
        cw = csv.writer(bfile, dialect='excel-tab')
        cw.writerow(header)
        for row in rows:
            cw.writerow(row)

    experiments = [ experiment.getExperimentById(exp_id) for exp_id in experiment_list ]
    downloadutils.experiment_metadata_to_tsv(experiments, metadata_filepath)

    downloadutils.zip_package([data_filepath, metadata_filepath], zip_filepath)

    exp = experiment.getExperimentById(exp_id, secure=False, check_ready=False)
    exp.delete()

    return success, errors

def create_temp_experiment(user_id, job_id):
    exp = experiment.Experiment()
    exp.name = 'temp experiment %d' % (random.randint(0,100000))
    exp.author = ''
    exp.description = ''
    exp.contact = ''
    exp.PMID=0
    exp.URL=''
    exp.published=0
    exp.ambiguity=0
    exp.experiment_id=None
    exp.dataset=''
    exp.volume=0
    exp.page_start=''
    exp.page_end=''
    exp.journal=''
    exp.publication_year=0
    exp.publication_month=''
    exp.public = 0
    exp.job_id = job_id
    exp.submitted_id = user_id
    exp.type='dataset'

    exp.saveExperiment()
    return exp.id

@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.dynamic_transaction_task
def batch_annotate_proteins(accessions, batch_id, user_id, job_id):
    notify_tasks.set_job_status.apply_async((job_id, 'started'))
    notify_tasks.set_job_stage.apply_async((job_id, 'initializing', 0))

    accession_dict = {}
    line_mapping = {}
    for i, acc in enumerate(accessions):
        accession_dict[acc] = set([i+1])
        line_mapping[i+1] = (acc, '')

    exp_id = create_temp_experiment(user_id, job_id)

    get_proteins_task = protein_tasks.get_proteins_from_external_databases.s(accession_dict, line_mapping, exp_id, job_id)
    get_protein_metadata_task = protein_tasks.query_protein_metadata.s(accession_dict, line_mapping, exp_id, job_id)
    annotate_proteins_task = annotate_proteins.s(accessions, batch_id, exp_id, user_id, job_id)
    notify_task = notify_tasks.finalize_batch_annotate_job.s(job_id)

    load_task = ( get_proteins_task | get_protein_metadata_task | annotate_proteins_task | notify_task )

    return load_task, (None,), None
