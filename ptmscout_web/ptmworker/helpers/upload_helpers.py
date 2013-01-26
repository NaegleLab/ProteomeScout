from celery.canvas import group
from ptmscout.database import protein, taxonomies, modifications, experiment, gene_expression, mutations, uniprot
from ptmscout.utils import uploadutils
from ptmscout.config import strings, settings
from ptmscout.database.modifications import NoSuchPeptide
from ptmworker.helpers import scansite_tools, quickgo_tools
import logging
import pickle
import re
import sys, os
import transaction
import traceback

log = logging.getLogger('ptmscout')

def transaction_task(fn):
    def ttask(*args):
        from ptmscout.database import DBSession
        log.debug("Running task: %s", fn.__name__)
        try:
            result = fn(*args)
            transaction.commit()
            return result
        except Exception:
            transaction.abort()
            raise
    ttask.__name__ = fn.__name__
    return ttask


def dynamic_transaction_task(fn):
    def ttask(*args):
        log.debug("Running task: %s", fn.__name__)
        try:
            result = fn(*args)
            transaction.commit()

            if result != None and len(result) == 3:
                new_task, task_arg, errback = result
                new_task.apply_async((task_arg,), link_error=errback )
        except Exception:
            transaction.abort()
            raise
    ttask.__name__ = fn.__name__
    return ttask


def find_activation_loops(prot):
    kinase_domain_names = set([ 'pkinase', 'pkinase_tyr' ])
    kinase_domains = [ d for d in prot.domains if d.label.lower() in kinase_domain_names ]
    cutoff_loop_size = 35

    start_motif = r"D[FPLY]G"
    stop_motif = r"[ASP][PILW][ED]"

    for d in kinase_domains:
        domain_seq = prot.sequence[d.start-1: d.stop]
        m1 = re.search(start_motif, domain_seq)

        if m1 == None:
            continue

        domain_seq = domain_seq[m1.end():]
        m2 = re.search(stop_motif, domain_seq)

        if m2 == None:
            continue

        domain_seq = domain_seq[:m2.start()]

        loop_start = d.start + m1.end()
        loop_end = d.start + m1.end() + m2.start() - 1

        label = strings.kinase_loop_name if len(domain_seq) <= 35 else strings.possible_kinase_name
        source = 'predicted'
        region = protein.ProteinRegion(label, source, loop_start, loop_end)

        prot.regions.append(region)



def check_ambiguity(measured_pep, species_name):
    for swissprot in uniprot.findPeptide(measured_pep.peptide, species_name):
        amb = modifications.PeptideAmbiguity(swissprot.accession, measured_pep.id)
        measured_pep.ambiguities.append(amb)


def parse_variants(acc, prot_seq, variants):
    new_mutations = []
    j = 0
    for mutantDict in variants:
        # for now we're only working with single mutants, but could expand
        # to double mutants in the future...
        if(mutantDict['type'] == "Substitution (single)"):
            new_mutation = mutations.Mutation(mutantDict['type'],
                    mutantDict['location'], mutantDict['original'],
                    mutantDict['mutant'], acc, mutantDict['notes'], None)

            if not new_mutation.consistent(prot_seq):
                log.info( "Loaded mutation does not match protein sequence (%d %s) %s -> %s" % (new_mutation.location, prot_seq[new_mutation.location-1], new_mutation.original, new_mutation.mutant) )
            else:
                new_mutations.append(new_mutation)

    return new_mutations

def get_go_annotation(goId, complete_go_terms, missing_terms):
    go_term = protein.getGoAnnotationById(goId)
    entry = None

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
        go_term.save()

    return go_term, entry

def query_missing_GO_terms(missing_terms, complete_go_terms):
    processed_missing = set()
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
        processed_missing.add(goId)

        for parent_goId in entry.is_a:
            if parent_goId not in complete_go_terms and \
                    protein.getGoAnnotationById(parent_goId) == None:
                missing_terms.append(parent_goId)


def report_errors(exp_id, errors, line_mapping):
    for e in errors:
        accession, peptide = line_mapping[e.row]
        experiment.createExperimentError(exp_id, e.row, accession, peptide, e.msg)


def get_strain_or_isolate(species):
    m = re.match(r"^(.*) \(((?:strain|isolate) .*)\)$", species)

    if m:
        species_root = m.group(1)
        strain = m.group(2)
        return species_root, strain

    return species, None

def get_taxonomic_lineage(species):
    species, strain = get_strain_or_isolate(species)
    taxon = taxonomies.getTaxonByName(species, strain)

    if taxon == None:
        return []
    
    taxonomic_lineage = []
    while taxon.parent != None:
        taxon = taxon.parent
        taxonomic_lineage.append(taxon.name)

    taxonomic_lineage.reverse()

    return taxonomic_lineage

def find_or_create_species(species):
    sp = taxonomies.getSpeciesByName(species)
    
    if(sp == None):
        species_root, strain = get_strain_or_isolate(species)
        tx = taxonomies.getTaxonByName(species_root, strain=strain)
        
        if tx == None:
            raise uploadutils.ParseError(None, None, "Species: " + species + " does not match any taxon node")
        
        sp = taxonomies.Species(species)
        sp.taxon_id = tx.node_id
        
    return sp


def create_accession_for_protein(prot, other_accessions):
    added_accessions = []

    for db, acc in other_accessions:
        db = db.lower()
        if not prot.hasAccession(acc):
            dbacc = protein.ProteinAccession()
            dbacc.type = db
            dbacc.value = acc
            prot.accessions.append(dbacc)
            added_accessions.append((db,acc))

    return added_accessions


def map_expression_probesets(prot):
    search_accessions = [ acc.value for acc in prot.accessions ]
    if prot.acc_gene != '' and prot.acc_gene != None:
        search_accessions.append(prot.acc_gene)
    
    probesets = gene_expression.getExpressionProbeSetsForProtein(search_accessions, prot.species_id)
    
    prot.expression_probes = []
    prot.expression_probes.extend(probesets)
    
    log.info("Loaded %d probesets for protein %s | %s", len(probesets), prot.accessions[0].value, str(prot.acc_gene))


def create_new_protein(name, gene, seq, species):
    prot = protein.Protein()
    prot.acc_gene = gene
    prot.name = name
    prot.sequence = seq
    prot.species = find_or_create_species(species)
    prot.species_id = prot.species.id
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
    pep_seq = pep_seq.strip()
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
        
        acc = row[acc_col.column_number].strip()
        pep = row[pep_col.column_number].strip()
        mods = row[mod_col.column_number].strip()
        
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
        
        mod_set = mod_map.get((acc,pep), set())
        mod_set.add(mods)
        mod_map[(acc,pep)] = mod_set

        run_data = data_runs.get((acc, pep, mods), {})
        
        series = []
        for d in data_cols:
            series.append(row[d.column_number].strip())
        for s in stddev_cols:
            series.append(row[s.column_number].strip())
            
        if run_col != None:
            run_data[ row[run_col.column_number].strip() ] = (line, series)
        else:
            run_data[ 'average' ] = (line, series)
        
        data_runs[(acc, pep, mods)] = run_data
        
    
    return accessions, peptides, mod_map, data_runs, errors, line_mapping
    

def extract_uniprot_accessions(accessions):
    uniprot_accs = []
    other_accs = []
    for acc in accessions:
        if(re.search('^[A-NR-Z]\d[A-Z]..\d([\.\-]\d+)?$', acc) != None):
            uniprot_accs.append(acc)
        elif(re.search('^[OPQ]\d...\d([\.\-]\d+)?$', acc) != None):
            uniprot_accs.append(acc)
        else:
            other_accs.append(acc)
    return uniprot_accs, other_accs

def group_critera(group, arg):
    if len(group) == 0:
        return True
    return group[-1][:6] == arg[:6]


def create_chunked_tasks_preserve_groups(task_args, MAX_BATCH_SIZE):
    tasks = []
    args = []

    small_groups = [[]]
    for arg in task_args:
        if group_critera(small_groups[-1], arg):
            small_groups[-1].append(arg)
        else:
            small_groups.append([arg])

    for g in small_groups:
        if len(args) + len(g) <= MAX_BATCH_SIZE:
            args = args + g
        else:
            tasks.append( args )
            args = g

    if len(args) > 0:
        tasks.append( args )

    return tasks


def create_chunked_tasks(task_args, MAX_BATCH_SIZE):
    tasks = []
    args = []
    
    for arg in task_args:
        args.append(arg)
        if len(args) == MAX_BATCH_SIZE:
            tasks.append( args )
            args = []
    if len(args) > 0:
        tasks.append( args )
        
    return tasks
    
def store_stage_input(exp_id, stage, result):
    stage = stage.replace(" ", "_")
    result_path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'e%d') % (exp_id)
    if not os.path.exists(result_path):
        os.mkdir(result_path)
    result_file = '%s.input' % (stage)

    writable = open(os.path.join(result_path, result_file), 'w')
    writable.write(pickle.dumps(result))
    writable.close()

def get_stage_input(exp_id, stage):
    stage = stage.replace(" ", "_")
    result_path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'e%d') % (exp_id)
    result_file = '%s.input' % (stage)

    file_path = os.path.join(result_path, result_file)
    if not os.path.exists(file_path):
        return None

    readable = open(file_path, 'r')
    result = pickle.loads(readable.read())
    readable.close()

    return result
