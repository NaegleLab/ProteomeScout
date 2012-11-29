from ptmscout.database import protein, taxonomies, modifications, experiment,\
    gene_expression
from ptmscout.utils import uploadutils
from ptmscout.config import strings
import logging
import sys
from ptmscout.database.modifications import NoSuchPeptide
from ptmworker import scansite_tools
import transaction
import traceback
from celery.canvas import group

log = logging.getLogger('ptmscout')


def logged_task(fn):
    def ttask(*args):
        log.debug("Running task: %s", fn.__name__)
        try:
            return fn(*args)
        except Exception:
            log.warning(traceback.format_exc())
    ttask.__name__ = fn.__name__
    return ttask


def transaction_task(fn):
    def ttask(*args):
        log.debug("Running task: %s", fn.__name__)
        try:
            result = fn(*args)
            transaction.commit()
            return result
        except Exception:
            log.warning(traceback.format_exc())
            transaction.abort()
    ttask.__name__ = fn.__name__
    return ttask


def dynamic_transaction_task(fn):
    def ttask(*args):
        log.debug("Running task: %s", fn.__name__)
        try:
            result = fn(*args)
            transaction.commit()
            
            if result != None and len(result) > 0:
                new_tasks = group(result)
                new_tasks.apply_async()
        except Exception:
            log.warning(traceback.format_exc())
            transaction.abort()
    ttask.__name__ = fn.__name__
    return ttask



def create_accession_for_protein(prot, other_accessions):
    for db, acc, _ in other_accessions:
        if not prot.hasAccession(acc):
            dbacc = protein.ProteinAccession()
            dbacc.type = db
            dbacc.value = acc
            prot.accessions.append(dbacc)


def map_expression_probesets(prot):
    search_accessions = [ acc.value for acc in prot.accessions ]
    if prot.acc_gene != '' and prot.acc_gene != None:
        search_accessions.append(prot.acc_gene)
    
    probesets = gene_expression.getExpressionProbeSetsForProtein(search_accessions, prot.species_id)
    
    prot.expression_probes.extend(probesets)
    
    log.info("Loaded %d probesets for protein %d | %s", len(probesets), prot.id, str(prot.acc_gene))


def report_errors(exp_id, errors, line_mapping):
    for e in errors:
        accession, peptide = line_mapping[e.row]
        experiment.createExperimentError(exp_id, e.row, accession, peptide, e.msg)


def mark_experiment(exp_id, status):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.status = status
    exp.saveExperiment()
    return exp

def find_or_create_species(species):
    sp = taxonomies.getSpeciesByName(species)
    
    if(sp == None):
        tx = taxonomies.getTaxonByName(species)
        
        if tx == None:
            raise uploadutils.ParseError(None, None, "Species: " + species + " does not match any taxon node")
        
        sp = taxonomies.Species()
        sp.name = species
        sp.taxon_id = tx.id
        
    return sp


def create_domains_for_protein(prot, domains, source, params):
    for domain in domains:
        dbdomain = protein.ProteinDomain()
        dbdomain.p_value = domain.p_value
        dbdomain.start = domain.start
        dbdomain.stop = domain.stop
        dbdomain.source = source
        dbdomain.version = domain.release
        dbdomain.label = domain.label
        dbdomain.params = params
        prot.domains.append(dbdomain)


def create_new_protein(name, gene, seq, species, accessions):
    log.info("Creating protein: %s SEQ: %s", str(accessions), seq)
    prot = protein.Protein()
    prot.acc_gene = gene
    prot.name = name
    prot.sequence = seq
    prot.species = find_or_create_species(species)
    
    for acc_type, acc in accessions:
        acc_db = protein.ProteinAccession()
        acc_db.type = acc_type
        acc_db.value = acc
        prot.accessions.append(acc_db)

    return prot

def get_related_proteins(prot_accessions, species):
    related_proteins = []
    for p in protein.getProteinsByAccession(prot_accessions, species):
        related_proteins.append(p)
    return related_proteins

def get_aligned_peptide_sequences(mod_sites, index, pep_seq, prot_seq):
    upper_case = pep_seq.upper()
    aligned_peptides = []
    
    for i in mod_sites:
        pep_site = i + index
        
        low_bound = max([pep_site-7, 0])
        high_bound = min([len(prot_seq), pep_site+8])

#        pep_tryps   = upper_case[:i] + pep_seq[i] + upper_case[i+1:]
        pep_aligned = prot_seq[low_bound:pep_site] + pep_seq[i] + prot_seq[pep_site+1:high_bound]
        
        if pep_site-7 < 0:
            pep_aligned = (" " * (7 - pep_site)) + pep_aligned
        if pep_site+8 > len(prot_seq):
            pep_aligned = pep_aligned + (" " * (pep_site + 8 - len(prot_seq)))
        
        pep_type = upper_case[i]
        
        aligned_peptides.append((pep_site+1, pep_aligned, pep_type))
    
    return aligned_peptides
    

def check_peptide_matches_protein_sequence(prot_seq, pep_seq):   
    index = prot_seq.find(pep_seq.upper())
    
    if index == -1:
        raise uploadutils.ParseError(None, None, strings.experiment_upload_warning_peptide_not_found_in_protein_sequence)
    
    return index

def parse_modifications(prot_seq, pep_seq, mods, taxonomy):
    index = check_peptide_matches_protein_sequence(prot_seq, pep_seq)
    mod_indices, mod_types = uploadutils.check_modification_type_matches_peptide(None, pep_seq, mods, taxonomy)
    aligned_sequences = get_aligned_peptide_sequences(mod_indices, index, pep_seq, prot_seq)
    
    return mod_types, aligned_sequences


def query_peptide_predictions(pep_seq, motif_class):
    log.info("Loading scansite predictions...")
    scansite_predictions = scansite_tools.get_scansite_motif(pep_seq, motif_class)
    
    db_predictions = []
    for scansite in scansite_predictions:
        pred = modifications.ScansitePrediction()
        pred.score = scansite.score
        pred.value = scansite.nickname
        pred.source = scansite.parse_source()
        db_predictions.append(pred)
#        log.info("%f %s %s", scansite.score, scansite.nickname, scansite.parse_source())
        
    return db_predictions

def get_peptide(prot_id, pep_site, peptide_sequence):
    upper_pep = peptide_sequence.upper()
    pep_type = upper_pep[7]
    
    created = False
    try:
        pep = modifications.getPeptideBySite(pep_site, pep_type, prot_id)
    except NoSuchPeptide:
        pep = modifications.Peptide()
        pep.pep_aligned = peptide_sequence
        pep.site_pos = pep_site
        pep.site_type = pep_type
        pep.protein_id = prot_id
        pep.save()
        created = True
    
    return pep, created
    

def insert_run_data(MSpeptide, line, units, series_header, run_name, series):
    for i in xrange(0, len(series_header)):
        try:
            tp, x = series_header[i]
            y = float(series[i])

            data = MSpeptide.getDataElement(run_name, tp, x)

            if data == None:
                data = experiment.ExperimentData()
                MSpeptide.data.append(data)

            data.run = run_name
            data.priority = i + 1
            data.type = tp
            data.units = units
            data.label = x
            data.value = y
        except Exception, e:
            log.warning("Error inserting data element on line %d: '%s' exc: %s", line, series[i], str(e))


def get_series_headers(session):
    headers = []
    for col in session.getColumns('data'):
        headers.append(('data', col.label))
    
    for col in session.getColumns('stddev'):
        headers.append(('stddev', col.label))
    
    return headers



def parse_datafile(session):
    accessions = {}
    peptides = {}
    mod_map = {}
    data_runs = {}
    line_mapping = {}
    errors = []
    
    acc_col = session.getColumns('accession')[0]
    pep_col = session.getColumns('peptide')[0]
    mod_col = session.getColumns('modification')[0]
    run_col = None
    
    found_cols = session.getColumns('run')
    if found_cols != []:
        run_col = found_cols[0]
    
    data_cols = session.getColumns('data')
    stddev_cols = session.getColumns('stddev')
    
    _, rows = uploadutils.load_header_and_data_rows(session.data_file, sys.maxint)
    
    keys = set([])

    line=0
    for row in rows:
        line+=1
        line_errors = uploadutils.check_data_row(line, row, acc_col, pep_col, mod_col, run_col, data_cols, stddev_cols, keys)
        
        acc = row[acc_col.column_number]
        pep = row[pep_col.column_number]
        mods = row[mod_col.column_number]
        
        line_mapping[line] = (acc, pep)
        
        if len(line_errors) > 0:
            errors.extend(line_errors)
            continue
        
        line_set = accessions.get(acc, [])
        line_set.append(line)
        accessions[acc] = line_set
        
        pep_set = peptides.get(acc, set())
        pep_set.add(pep)
        peptides[acc] = pep_set
        
        mod_map[(acc,pep)] = mods
        
        run_data = data_runs.get((acc, pep), {})
        
        series = []
        for d in data_cols:
            series.append(row[d.column_number])
        for s in stddev_cols:
            series.append(row[s.column_number])
            
        if run_col != None:
            run_data[ row[run_col.column_number] ] = (line, series)
        else:
            run_data[ 'average' ] = (line, series)
        
        data_runs[(acc,pep)] = run_data
        
    
    return accessions, peptides, mod_map, data_runs, errors, line_mapping
    

def create_chunked_tasks(task_method, task_args, MAX_BATCH_SIZE):
    tasks = []
    args = []
    
    for arg in task_args:
        args.append(arg)
        if len(args) == MAX_BATCH_SIZE:
            tasks.append( task_method.s(args) )
            args = []
    if len(args) > 0:
        tasks.append( task_method.s(args) )
        
    return tasks
    

