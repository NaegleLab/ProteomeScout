from ptmscout.database import protein, taxonomies, modifications, experiment
from ptmscout.utils import protein_utils, uploadutils
from ptmscout.config import strings
import logging
import sys
from ptmscout.database.modifications import NoSuchPeptide
from ptmworker import scansite_tools

log = logging.getLogger('ptmscout')

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
    log.debug("Creating protein: %s SEQ: %s", str(accessions), seq)
    prot = protein.Protein()
    prot.acc_gene = gene
    prot.name = name
    prot.sequence = seq
    prot.species = find_or_create_species(species)
    
    for prot_acc in accessions:
        acc_db = protein.ProteinAccession()
        acc_db.type = protein_utils.get_accession_type(prot_acc)
        acc_db.value = prot_acc
        prot.accessions.append(acc_db)

    return prot


def find_protein(seq, prot_accessions, species):
    for p in protein.getProteinsByAccession(prot_accessions, species):
        if p.sequence.lower() == seq.lower():
            return p

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
    
    mod_map = {}
    
    for i in xrange(0, len(mod_indices)):
        mod_map[mod_indices[i]] = (mod_types[i], aligned_sequences[i])
    
    return mod_map


def query_peptide_predictions(pep_seq, motif_class):
    scansite_predictions = scansite_tools.get_scansite_motif(pep_seq, motif_class)
    
    db_predictions = []
    for scansite in scansite_predictions:
        pred = modifications.ScansitePrediction()
        pred.score = scansite.score
        pred.value = scansite.nickname
        pred.source = scansite.parse_source()
        db_predictions.append(pred)
        
    return db_predictions

def get_peptide(prot_id, pep_site, peptide_sequence):
    upper_pep = peptide_sequence.upper()
    pep_type = upper_pep[8]
    
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
            data = experiment.ExperimentData()
            
            tp, x = series_header[i]
            
            y = float(series[i])
            
            data.run = run_name
            data.priority = i + 1
            data.type = tp
            data.units = units
            data.label = x
            data.value = y
            MSpeptide.data.append(data)
        except Exception, e:
            log.debug("Error inserting data element on line %d: '%s' exc: %s", line, series[i], str(e))


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
    
