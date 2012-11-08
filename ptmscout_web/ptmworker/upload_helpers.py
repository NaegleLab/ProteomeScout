from ptmscout.database import protein, taxonomies, modifications, experiment
from ptmscout.utils import protein_utils, uploadutils
from ptmscout.config import settings, strings
import logging
from geeneus import Proteome
import sys

log = logging.getLogger('ptmscout')


def mark_experiment(exp_id, status):
    exp = experiment.getExperimentById(exp_id, check_ready=False, secure=False)
    exp.status = status
    exp.saveExperiment()
    return exp

def get_protein_information(pm, acc):
    seq = pm.get_protein_sequence(acc)
    
    if seq == "":
        raise uploadutils.ParseError(None, None, "Unable to fetch protein accession: %s" % (acc))
    
    name = pm.get_protein_name(acc)
    gene = pm.get_geneID(acc)
    
    XML = pm.get_raw_xml(acc)

    taxonomy = XML[0]['GBSeq_taxonomy']
    taxonomy = set([ t.strip() for t in taxonomy.split(';') ])
    
    species = XML[0]['GBSeq_organism']
    
    prot_accessions = [ XML[0]['GBSeq_primary-accession'] ]
    prot_accessions.extend( XML[0]['GBSeq_other-seqids'] )
    
    return name, gene, taxonomy, species, prot_accessions, seq


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

def create_new_protein(name, gene, seq, species, accessions):
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

def find_protein(name, gene, seq, prot_accessions, species):
    for p in protein.getProteinsByAccession(prot_accessions, species):
        if p.sequence.lower() == seq.lower():
            return p
    
    log.debug("Creating protein: %s SEQ: %s", str(prot_accessions), seq)
    return create_new_protein(name, gene, seq, species, prot_accessions)


def load_proteins(accessions, pm, prot_map):
    log.debug("Querying for proteins: %s", str(accessions))
    pm.batch_get_protein_sequence(accessions)


def get_proteins_from_ncbi(accessions, MAX_BATCH_SIZE):
    pm = Proteome.ProteinManager(settings.adminEmail)

    query_job_args = [[]]
    for acc in accessions:
        query_job_args[-1].append(acc)
        
        if len(query_job_args[-1]) == MAX_BATCH_SIZE:
            query_job_args.append([])

        
    prot_map = {}
    for qjob in query_job_args:
        load_proteins(qjob, pm, prot_map)
        

    errors = []
    for acc in accessions:
        try:
            prot_map[acc] = get_protein_information(pm, acc)
        except uploadutils.ParseError, e:
            for line in accessions[acc]:
                e.row = line
                errors.append(e)
    
            
    return prot_map, errors


def get_aligned_peptide_sequences(mod_sites, index, pep_seq, prot_seq):
    upper_case = pep_seq.upper()
    
    aligned_peptides = []
    
    for i in mod_sites:
        pep_site = i + index
        
#        pep_tryps   = upper_case[:i] + pep_seq[i] + upper_case[i+1:]
        pep_aligned = prot_seq[pep_site-7:pep_site] + pep_seq[i] + prot_seq[pep_site+1:pep_site+8]
        pep_type = upper_case[i]
        
        aligned_peptides.append((pep_site, pep_aligned, pep_type))
    
    return aligned_peptides
    

def create_modifications(protein_id, prot_seq, pep_seq, mods, taxonomy):
    created_mods = []
    
    prot_seq = prot_seq.upper()
    
    index = prot_seq.find(pep_seq.upper())
    if index == -1:
        raise uploadutils.ParseError(None, None, strings.experiment_upload_warning_peptide_not_found_in_protein_sequence)
    
    mod_indices, mod_types = uploadutils.check_modification_type_matches_peptide(None, pep_seq, mods, taxonomy)
    aligned_sequences = get_aligned_peptide_sequences(mod_indices, index, pep_seq, prot_seq)
    
    for i in xrange(0, len(mod_indices)):
        mod = mod_types[i]
        pep_site, pep_aligned, pep_type = aligned_sequences[i]
        
        try:
            pep_mod = modifications.getModificationBySite(pep_site, pep_type, protein_id, mod.id)
        except modifications.NoSuchModification:
            pep_mod = modifications.Phosphopep()
            pep_mod.pep_aligned = pep_aligned
            pep_mod.pep_tryps = ""
            pep_mod.site_pos = pep_site
            pep_mod.site_type = pep_type
            pep_mod.protein_id = protein_id
            pep_mod.mod_id = mod.id
        
        created_mods.append(pep_mod)
        
    return created_mods
    

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
            data.MS = MSpeptide
            
            data.save()
        except Exception, e:
            log.debug("Error inserting data element on line %d: '%s' exc: %s", line, series[i], str(e))



def parse_datafile(session):
    accessions = {}
    peptides = {}
    mod_map = {}
    data_runs = {}
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
    
    _, rows = uploadutils.load_header_and_data_rows(session, sys.maxint)
    
    keys = set([])

    line = 1
    for row in rows:
        line_errors = uploadutils.check_data_row(line, row, acc_col, pep_col, mod_col, run_col, data_cols, stddev_cols, keys)
        
        if len(line_errors) > 0:
            errors.extend(line_errors)
            continue
        
        acc = row[acc_col.column_number]
        pep = row[pep_col.column_number]
        mods = row[mod_col.column_number]
        
        acc_lines = accessions.get(acc, [])
        acc_lines.append(line)
        accessions[acc] = acc_lines
        
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
        
        line+=1
    
    return accessions, peptides, mod_map, data_runs, errors
    

def get_series_headers(session):
    headers = []
    for col in session.getColumns('data'):
        headers.append(('data', col.label))
    
    for col in session.getColumns('stddev'):
        headers.append(('stddev', col.label))
    
    return headers
