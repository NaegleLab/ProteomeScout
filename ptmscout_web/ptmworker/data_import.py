import celery
import logging 
from ptmworker import upload_helpers, entrez_tools, pfam_tools, quickgo_tools
from ptmscout.database import upload, modifications, protein
from celery.canvas import group

log = logging.getLogger('ptmscout')

PFAM_DEFAULT_CUTOFF = 0.00001
MAX_NCBI_BATCH_SIZE = 1000
MAX_QUICKGO_BATCH_SIZE = 100

@celery.task
def finalize_import(exp_id, user_email, application_url):
    pass


@celery.task(rate_limit='3/s')
def load_new_peptide(prot_id, pep_id, pep_seq, pep_site, taxonomy):
    pep = modifications.getPeptideById(pep_id)
    pep.protein_domain = protein.getProteinDomain(prot_id, pep_site)

    motif_class = None
    if 'mammalia' in taxonomy:
        motif_class="MAMMALIAN"
    elif 'saccharomycotina' in taxonomy:
        motif_class="YEAST"
    elif 'saccharomyces' in taxonomy:
        motif_class="YEAST"

    if motif_class != None:
        pep.predictions = upload_helpers.query_peptide_predictions(pep_seq, motif_class)
        
    pep.save()
    
    

@celery.task
def load_peptide_modification(protein_info, exp_id, pep_seq, mods, units, series_header, runs):
    protein_id, protein_sequence, taxonomy = protein_info
    
    
    pep_measurement = modifications.MeasuredPeptide()
    
    pep_measurement.experiment_id = exp_id
    pep_measurement.peptide = pep_seq
    pep_measurement.protein_id = protein_id
        
    mod_map = upload_helpers.parse_modifications(protein_sequence, pep_seq, mods, taxonomy)
    
    new_peptide_tasks = []
    
    for site_pos in mod_map:
        mod_type, pep_aligned = mod_map[site_pos]
        
        pep, created = upload_helpers.get_peptide(protein_id, site_pos, pep_aligned)
        
        pepmod = modifications.PeptideModification()
        pepmod.modification = mod_type
        pepmod.peptide = pep
        pep_measurement.peptides.append(pepmod)
        
        for line, run_name, series in runs:
            upload_helpers.insert_run_data(pep_measurement, line, units, series_header, run_name, series)
        
        if created:
            new_peptide_tasks.append( load_new_peptide.s(protein_id, pep.id, pep_aligned, site_pos, taxonomy) )

    pep_measurement.save()
    
    if len(new_peptide_tasks) > 0:
        pep_tasks = group(new_peptide_tasks)
        pep_tasks.apply_async()


@celery.task
def load_protein(accession, protein_information):
    _a, _a, taxonomy, species, accessions, _d, seq = protein_information
    prot = upload_helpers.find_protein(seq, accessions, species)
    
    return prot.id, prot.sequence, taxonomy


@celery.task(rate_limit='3/s')
def load_new_protein(accession, protein_information):
    _a, _a, taxonomy, species, accessions, domains, seq = protein_information
    prot = upload_helpers.find_protein(seq, accessions, species)
    
    if len(domains) == 0:
        domains = pfam_tools.get_computed_pfam_domains(prot.sequence, PFAM_DEFAULT_CUTOFF)
        upload_helpers.create_domains_for_protein(prot, domains, "PARSED PFAM", "ALL")
    else:
        upload_helpers.create_domains_for_protein(prot, domains, "COMPUTED PFAM", "pval=%f" % (PFAM_DEFAULT_CUTOFF))
    
    # load additional protein accessions if available
    
    
    prot.save()
    
    return prot.id, prot.sequence, taxonomy


@celery.task
def create_missing_GO_annotations(GO_annotation_maps, protein_ids):
    aggregate_GO_map = {}
    for map in GO_annotation_maps:
        for acc in map:
            aggregate_GO_map[acc] = map[acc]
        
    

@celery.task(rate_limit='3/s')
def get_GO_annotations(protein_accessions):
    go_terms, _ = quickgo_tools.batch_get_GO_annotations(protein_accessions)
    return go_terms


def create_missing_proteins(prot_map, missing_proteins):
    protein_id_map = []
    
    for acc in missing_proteins:
        name, gene, _t, species, accessions, _d, seq = prot_map[acc]
        prot = upload_helpers.create_new_protein(name, gene, seq, species, accessions)
        prot.saveProtein()
        protein_id_map[acc] = prot.id
        
    return protein_id_map

def create_GO_import_tasks(protein_map, new_protein_ids):
    GO_annotation_tasks = upload_helpers.create_chunked_tasks(get_GO_annotations, protein_map, MAX_QUICKGO_BATCH_SIZE)
    GO_term_task = create_missing_GO_annotations.s(new_protein_ids)
    
    return ( group(GO_annotation_tasks) | GO_term_task )


def create_protein_import_tasks(prot_map, missing_proteins, parsed_datafile, headers):
    for acc in prot_map:
        pep_tasks = []
        
        
        
        



@celery.task(rate_limit='3/s')
def get_ncbi_proteins(protein_accessions):
    prot_map, errors = entrez_tools.get_proteins_from_ncbi(protein_accessions)
    return prot_map, errors



@celery.task
def launch_loader_tasks(ncbi_results, parsed_datafile, headers, session_info):
    #combine the results
    aggregate_errors = []
    aggregate_protein_map = {}
    for prot_map, errors in ncbi_results:
        aggregate_errors.extend(errors)
        for acc in prot_map:
            aggregate_protein_map[acc] = prot_map[acc]

    #list the missing proteins
    missing_proteins = set()
    for acc in aggregate_protein_map:
        _n, _g, _t, species, prot_accessions, _d, seq = aggregate_protein_map[acc]
        if upload_helpers.find_protein(seq, prot_accessions, species) == None:
            missing_proteins.add(acc)

    #create entries for the missing proteins
    new_protein_ids = create_missing_proteins(aggregate_protein_map, missing_proteins)
    
    #run go import for missing proteins
    GO_task = create_GO_import_tasks(aggregate_protein_map, new_protein_ids)
    
    protein_tasks = create_protein_import_tasks(aggregate_protein_map, missing_proteins, parsed_datafile, headers)
    protein_tasks.append(GO_task)
    
    exp_id, _s, user_email, application_url = session_info
    import_tasks = ( group(protein_tasks) | finalize_import(exp_id, user_email, application_url) )
    
    import_tasks.apply_async()


@celery.task
def start_import(exp_id, session_id, user_email, application_url):
    log.debug("starting import")
    upload_helpers.mark_experiment(exp_id, 'loading')
    
    session = upload.getSessionById(session_id, secure=False)
    session_info = (exp_id, session_id, user_email, application_url)
    
    accessions, peptides, mod_map, data_runs, errors, line_mapping = upload_helpers.parse_datafile(session)
    parsed_datafile = (accessions, peptides, mod_map, data_runs, errors, line_mapping)
    
    headers = upload_helpers.get_series_headers(session)
    
    ncbi_tasks = upload_helpers.create_chunked_tasks(get_ncbi_proteins, accessions, MAX_NCBI_BATCH_SIZE)
    
    
    load_task = ( group(ncbi_tasks) | launch_loader_tasks.s(parsed_datafile, headers, session_info) )
    load_task.apply_async()    