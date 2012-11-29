import celery
import logging 
from ptmworker import upload_helpers, entrez_tools, pfam_tools, quickgo_tools, picr_tools
from ptmscout.database import upload, modifications, protein, experiment
from celery.canvas import group
from ptmscout.config import strings, settings
from ptmscout.utils import mail, uploadutils

log = logging.getLogger('ptmscout')

PFAM_DEFAULT_CUTOFF = 0.00001
MAX_NCBI_BATCH_SIZE = 500
MAX_QUICKGO_BATCH_SIZE = 500

@celery.task
@upload_helpers.transaction_task
def finalize_experiment_error_state(exp_id):
    exp = upload_helpers.mark_experiment(exp_id, 'error')
    
    subject = strings.experiment_upload_failed_subject
    message = strings.experiment_upload_failed_message % (exp.name)
    
    mail.celery_send_mail([settings.adminEmail], subject, message)

@celery.task
@upload_helpers.transaction_task
def finalize_import(exp_id, user_email, application_url):
    exp = upload_helpers.mark_experiment(exp_id, 'loaded')

    peptides = modifications.getMeasuredPeptidesByExperiment(exp_id, secure=False, check_ready=False)
    proteins = set([ pep.protein_id for pep in peptides ])
    error_log_url = "%s/experiments/%d/errors" % (application_url, exp_id)
    
    subject = strings.experiment_upload_finished_subject
    message = strings.experiment_upload_finished_message % (exp.name, len(peptides), len(proteins), len(exp.errors), error_log_url)
    
    mail.celery_send_mail([user_email], subject, message)


@celery.task(rate_limit='5/s')
@upload_helpers.transaction_task
def load_new_peptide(prot_id, site_pos, pep_seq, taxonomy):
    pep, _ = upload_helpers.get_peptide(prot_id, site_pos, pep_seq)
    pep.protein_domain = protein.getProteinDomain(prot_id, site_pos)

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
@upload_helpers.dynamic_transaction_task
def load_peptide_modification(protein_info, exp_id, pep_seq, mods, units, series_header, runs):
    protein_id, protein_accession, protein_sequence, taxonomy = protein_info
    
    try:
        mod_types, aligned_sequences = upload_helpers.parse_modifications(protein_sequence, pep_seq, mods, taxonomy)

        pep_measurement = modifications.MeasuredPeptide()
        
        pep_measurement.experiment_id = exp_id
        pep_measurement.peptide = pep_seq
        pep_measurement.protein_id = protein_id
        
        new_peptide_tasks = []
        
        for i in xrange(0, len(aligned_sequences)):
            mod_type = mod_types[i] 
            site_pos, pep_sequence, _ = aligned_sequences[i]
            
            pep, created = upload_helpers.get_peptide(protein_id, site_pos, pep_sequence)
            
            pepmod = modifications.PeptideModification()
            pepmod.modification = mod_type
            pepmod.peptide = pep
            pep_measurement.peptides.append(pepmod)
            
            for line, run_name, series in runs:
                upload_helpers.insert_run_data(pep_measurement, line, units, series_header, run_name, series)
            
            if created:
                new_peptide_tasks.append( load_new_peptide.s(protein_id, site_pos, pep_sequence, taxonomy) )
    
        pep_measurement.save()
        
        return new_peptide_tasks

    except uploadutils.ParseError, pe:
        for line, _rn, _s in runs:
            experiment.createExperimentError(exp_id, line, protein_accession, pep_seq, pe.msg)


@celery.task
@upload_helpers.logged_task
def load_protein(accession, protein_information):
    _a, _a, taxonomy, species, _a, _d, seq = protein_information
    prot = protein.getProteinBySequence(seq, species)
    
    return prot.id, accession, prot.sequence, taxonomy


@celery.task(rate_limit='5/s')
@upload_helpers.transaction_task
def load_new_protein(accession, protein_information):
    _a, _a, taxonomy, species, _accessions, domains, seq = protein_information
    prot = protein.getProteinBySequence(seq, species)
    
    if len(domains) == 0:
        domains = pfam_tools.get_computed_pfam_domains(prot.sequence, PFAM_DEFAULT_CUTOFF)
        upload_helpers.create_domains_for_protein(prot, domains, "PARSED PFAM", "ALL")
    else:
        upload_helpers.create_domains_for_protein(prot, domains, "COMPUTED PFAM", "pval=%f" % (PFAM_DEFAULT_CUTOFF))
    
    # load additional protein accessions if available
    other_accessions = picr_tools.get_picr(accession, prot.species.taxon_id)
    upload_helpers.create_accession_for_protein(prot, other_accessions)
    
    upload_helpers.map_expression_probesets(prot)
    
    prot.saveProtein()
    
    return prot.id, accession, prot.sequence, taxonomy


def create_protein_import_tasks(prot_map, missing_proteins, parsed_datafile, headers, units, exp_id):
    _a, peptides, mod_map, data_runs, _l = parsed_datafile
    
    prot_tasks = []
    for acc in prot_map:
        pep_tasks = []
        
        for pep in peptides[acc]:
            key = (acc, pep)
            mod_str = mod_map[key]
            
            run_tasks = []
            for run_name in data_runs[key]:
                line, series = data_runs[key][run_name]
                run_tasks.append( (line, run_name, series) )
            
            pep_tasks.append( load_peptide_modification.s(exp_id, pep, mod_str, units, headers, run_tasks) )
        
        if acc in missing_proteins:
            prot_tasks.append( ( load_new_protein.s( acc, prot_map[acc] ) | group(pep_tasks) ) )
        else:
            prot_tasks.append( ( load_protein.s( acc, prot_map[acc] ) | group(pep_tasks) ) )
    
    return prot_tasks


@celery.task
@upload_helpers.transaction_task
def create_hierarchy_for_missing_GO_annotations(created_entries):
    links = 0
    for entry in created_entries:
        go_term = protein.getGoAnnotationById(entry.goId)
        
        for parent_goId in entry.is_a:
            parent_term = protein.getGoAnnotationById(parent_goId)
            
            if parent_term == None:
                log.warn("Parent term '%s' not found in GO annotations!", parent_goId)
            else:
                parent_term.children.append(go_term)
                parent_term.save()
                links += 1
    
    log.info("Created: %d GO entries with %d edges", len(created_entries), links)        


@celery.task
@upload_helpers.transaction_task
def create_missing_GO_annotations(GO_annotation_map, protein_ids):
    created_go_entries = []
    created, assigned = 0, 0
    
    complete_go_terms = set()
    for acc in GO_annotation_map:
        for goId in GO_annotation_map[acc]:
            complete_go_terms.add(goId)


    missing_terms = set()

    for acc in protein_ids:
        go_annotations = GO_annotation_map[acc]
        protein_id = protein_ids[acc]
        

        for goId in go_annotations:
            dateAdded = go_annotations[goId]
            go_term = protein.getGoAnnotationById(goId)
            
            if go_term == None:
                go_term = protein.GeneOntology()
                version, entry = quickgo_tools.get_GO_term(goId)

                for parent_goId in entry.is_a:
                    if parent_goId not in complete_go_terms and \
                            protein.getGoAnnotationById(parent_goId) == None:
                        missing_terms.add(parent_goId)

                go_term.GO = entry.goId
                go_term.term = entry.goName
                go_term.aspect = entry.goFunction
                go_term.version = version
                created_go_entries.append(entry)
                created+=1
                
            goe = protein.GeneOntologyEntry()
            goe.GO_term = go_term
            goe.protein_id = protein_id
            goe.date = dateAdded
            goe.save()
            assigned+=1
    
    processed_missing = set()
    missing_terms = list(missing_terms)
    while len(missing_terms) > 0:
        goId = missing_terms.pop(0)
        if goId in processed_missing:
            continue

        go_term = protein.GeneOntology()
        version, entry = quickgo_tools.get_GO_term(goId)

        go_term.GO = entry.goId
        go_term.term = entry.goName
        go_term.aspect = entry.goFunction
        go_term.version = version
        go_term.save()
        created_go_entries.append(entry)
        created+=1
        processed_missing.add(goId)

        for parent_goId in entry.is_a:
            if parent_goId not in complete_go_terms and \
                    protein.getGoAnnotationById(parent_goId) == None:
                missing_terms.append(parent_goId)


    log.info("Assigned %d terms, created %d terms", assigned, created)
    return created_go_entries


@celery.task
@upload_helpers.logged_task
def aggregate_GO_annotations(GO_annotation_maps):
    aggregate_GO_map = {}
    for GOmap in GO_annotation_maps:
        for acc in GOmap:
            aggregate_GO_map[acc] = GOmap[acc]
    
    return aggregate_GO_map

@celery.task(rate_limit='3/s')
@upload_helpers.logged_task
def get_GO_annotations(protein_accessions):
    go_terms, _ = quickgo_tools.batch_get_GO_annotations(protein_accessions)
    return go_terms


def create_GO_import_tasks(protein_map, new_protein_ids):
    GO_annotation_tasks = upload_helpers.create_chunked_tasks(get_GO_annotations, protein_map.keys(), MAX_QUICKGO_BATCH_SIZE)
    
    GO_aggregate_task = aggregate_GO_annotations.s()
    GO_term_task = create_missing_GO_annotations.s(new_protein_ids)
    GO_hierarchy_task = create_hierarchy_for_missing_GO_annotations.s()
    
    return ( group(GO_annotation_tasks) | GO_aggregate_task | GO_term_task | GO_hierarchy_task )


@celery.task(rate_limit='3/s')
@upload_helpers.logged_task
def get_ncbi_proteins(protein_accessions):
    log.info("Getting records for %d accessions", len(protein_accessions))
    prot_map, errors = entrez_tools.get_proteins_from_ncbi(protein_accessions)
    
    return prot_map, errors


@celery.task
@upload_helpers.logged_task
def aggregate_ncbi_results(ncbi_results, exp_id, accessions, line_mappings):
    
    #combine the results
    aggregate_errors = []
    aggregate_protein_map = {}
    for prot_map, errors in ncbi_results:
        aggregate_errors.extend(errors)
        for acc in prot_map:
            aggregate_protein_map[acc] = prot_map[acc]
    
    log.info("Detected %d errors", len(aggregate_errors))
    
    #report errors
    for error in aggregate_errors:
        for line in accessions[error.acc]:
            accession, peptide = line_mappings[line]
            experiment.createExperimentError(exp_id, line, accession, peptide, strings.experiment_upload_warning_accession_not_found)
    
    return aggregate_protein_map


@celery.task
@upload_helpers.transaction_task
def create_missing_proteins(protein_map):
    #list the missing proteins
    missing_proteins = set()
    for acc in protein_map:
        _n, _g, _t, species, _a, _d, seq = protein_map[acc]
        if protein.getProteinBySequence(seq, species) == None:
            missing_proteins.add(acc)

    #create entries for the missing proteins
    protein_id_map = {}
    for acc in missing_proteins:
        name, gene, _t, species, accessions, _d, seq = protein_map[acc]
        prot = upload_helpers.create_new_protein(name, gene, seq, species, accessions)
        prot.saveProtein()
        protein_id_map[acc] = prot.id
    
    return protein_map, missing_proteins, protein_id_map

@celery.task
@upload_helpers.logged_task
def launch_loader_tasks(ncbi_result, parsed_datafile, headers, session_info):
    log.info("Launching loader...")
    protein_map, missing_proteins, new_protein_ids = ncbi_result
    exp_id, _s, user_email, application_url, units = session_info
    
    #run go import for missing proteins
    GO_task = create_GO_import_tasks(protein_map, new_protein_ids)
    
    protein_tasks = create_protein_import_tasks(protein_map, missing_proteins, parsed_datafile, headers, units, exp_id)
    protein_tasks.append(GO_task)
    
    import_tasks = ( group(protein_tasks) | finalize_import.si(exp_id, user_email, application_url) )
    
    import_tasks.apply_async(link_error=finalize_experiment_error_state.s(exp_id))


@celery.task
@upload_helpers.transaction_task
def start_import(exp_id, session_id, user_email, application_url):
    exp = upload_helpers.mark_experiment(exp_id, 'loading')
    
    log.info("Loading session info...")
    session = upload.getSessionById(session_id, secure=False)
    session_info = (exp_id, session_id, user_email, application_url, session.units)
    
    log.info("Loading data file...")
    accessions, peptides, mod_map, data_runs, errors, line_mapping = upload_helpers.parse_datafile(session)

    exp.num_measured_peptides = len(mod_map)
    exp.saveExperiment()

    log.info("Reporting data file errors...")
    upload_helpers.report_errors(exp_id, errors, line_mapping)
    
    if len(accessions) == 0:
        log.info("Nothing to do: all proteins were rejected!")
        finalize_import.apply_async((exp_id, user_email, application_url))
    else:
        headers = upload_helpers.get_series_headers(session)
        parsed_datafile = (accessions, peptides, mod_map, data_runs, line_mapping)
        
        log.info("Running tasks...")
        ncbi_tasks = upload_helpers.create_chunked_tasks(get_ncbi_proteins, accessions.keys(), MAX_NCBI_BATCH_SIZE)
        
        load_task = ( group(ncbi_tasks) | aggregate_ncbi_results.s(exp_id, accessions, line_mapping) | create_missing_proteins.s() | launch_loader_tasks.s(parsed_datafile, headers, session_info) )
        load_task.apply_async(link_error=finalize_experiment_error_state.s(exp_id))
    
        log.info("Tasks started... now we wait")
    