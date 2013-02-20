import celery
import logging
from ptmworker import notify_tasks
from ptmworker.helpers import upload_helpers, quickgo_tools
from ptmscout.database import protein

log = logging.getLogger('ptmscout')

MAX_QUICKGO_BATCH_SIZE = 500

def create_hierarchies(created_entries):
    links = 0
    for goId in created_entries:
        entry = created_entries[goId]
        go_term = protein.getGoAnnotationById(goId)
        
        for parent_goId in entry.is_a:
            parent_term = protein.getGoAnnotationById(parent_goId)
            
            if parent_term == None:
                log.warn("Parent term '%s' not found in GO annotations!", parent_goId)
            else:
                parent_term.children.append(go_term)
                parent_term.save()
                links += 1

    log.info("Created: %d new GO is_a relationships", links)

def get_go_annotation(goId, created_entries):
    go_term = protein.getGoAnnotationById(goId)
    entry = None

    if go_term == None:
        go_term = protein.GeneOntology()
        version, entry = quickgo_tools.get_GO_term(goId)

        created_entries[goId] = entry

        go_term.GO = entry.goId
        go_term.term = entry.goName
        go_term.aspect = entry.goFunction
        go_term.version = version
        go_term.save()

    return go_term

def query_missing_GO_terms(created_entries):
    missing_terms = set()

    for goId in created_entries:
        entry = created_entries[goId]
        for parent_goId in entry.is_a:
            if protein.getGoAnnotationById(parent_goId) == None:
                missing_terms.add(parent_goId)

    while len(missing_terms) > 0:
        goId = missing_terms.pop()
        if goId in created_entries:
            continue

        go_term = protein.GeneOntology()
        version, entry = quickgo_tools.get_GO_term(goId)

        go_term.GO = entry.goId
        go_term.term = entry.goName
        go_term.aspect = entry.goFunction
        go_term.version = version
        go_term.save()
        created_entries[goId] = entry

        for parent_goId in entry.is_a:
            if protein.getGoAnnotationById(parent_goId) == None:
                missing_terms.add(parent_goId)

def get_primary_accessions(prot):
    uniprot_accessions = []
    refseq_accessions = []
    genbank_accessions = []

    for acc in prot.accessions:
        tp = acc.type.lower()
        if tp == 'uniprot' or tp == 'swissprot':
            uniprot_accessions.append(acc.value)
        if tp == 'refseq':
            refseq_accessions.append(acc.value)
        if tp == 'genbank':
            genbank_accessions.append(acc.value)

    primary_accessions = sorted(uniprot_accessions, key=lambda item: len(item))
    if len(primary_accessions) > 0:
        return primary_accessions
    if len(refseq_accessions) > 0:
        return refseq_accessions
    return genbank_accessions


def annotate_proteins(acc, annotations, protein_ids, protein_map, created_entries):
    assigned = 0
    duplicates = 0
    terms = set()
    for goId in annotations:
        terms.add(goId)

        go_term = get_go_annotation(goId, created_entries)
        dateAdded = annotations[goId]

        unassigned_proteins = [ protein_map[pid] for pid in protein_ids if not protein_map[pid].hasGoTerm(goId) ]
        assigned_proteins = [ protein_map[pid] for pid in protein_ids if protein_map[pid].hasGoTerm(goId) ]
        for prot in unassigned_proteins:
            prot.addGoTerm(go_term, dateAdded)
            assigned += 1
        duplicates += len(assigned_proteins)

    log.info( "Assigned %d new annotations to accession '%s' (%d proteins) from %d terms. Did not add %d duplicate annotations.", assigned, acc, len(protein_ids), len(terms), duplicates )

def assign_annotations(quickgo_result, protein_accessions, protein_map):
    created_entries = {}
    protein_ids = {}
    for pid in protein_accessions:
        for acc in protein_accessions[pid]:
            pids = protein_ids.get(acc, set())
            pids.add(pid)
            protein_ids[acc] = pids

    for acc in quickgo_result:
        if acc not in protein_ids:
            log.warning("Accession %s not found in stored protein accessions", acc)
        elif len(quickgo_result[acc]) > 0:
            annotate_proteins(acc, quickgo_result[acc], protein_ids[acc], protein_map, created_entries)
        else:
            log.info( "No annotations for accession '%s'", acc )

    for pid in protein_map:
        protein_map[pid].saveProtein()

    return created_entries

def build_accession_map(protein_ids):
    protein_accession_map = {}
    protein_map = {}

    for protein_id in protein_ids:
        prot = protein.getProteinById(protein_id)
        protein_map[protein_id] = prot
        protein_accession_map[protein_id] = get_primary_accessions(prot)

    return protein_accession_map, protein_map


@celery.task
@upload_helpers.transaction_task
def import_go_terms(protein_result, exp_id):
    protein_map, new_protein_ids = protein_result

    protein_accessions, protein_id_map = build_accession_map(new_protein_ids.values())
    query_accessions = set([ acc for pid in protein_accessions for acc in protein_accessions[pid] ])

    GO_map = {}
    queries = upload_helpers.create_chunked_tasks(query_accessions, MAX_QUICKGO_BATCH_SIZE)
    max_progress = len(queries) + 1

    upload_helpers.store_stage_input(exp_id, 'GO terms', protein_result)
    notify_tasks.set_loading_stage.apply_async((exp_id, 'GO terms', max_progress))

    i=0
    for qaccessions in queries:
        go_terms, _ = quickgo_tools.batch_get_GO_annotations(qaccessions)
        GO_map.update( go_terms )
        i+=1
        notify_tasks.set_progress.apply_async((exp_id, i, max_progress))

    log.info("Unexpected accessions returned from QuickGO: %s", str(set(GO_map.keys()) - query_accessions))

    created_entries = assign_annotations(GO_map, protein_accessions, protein_id_map)

    log.info( "Created a total of %d new GO terms", len(created_entries) )

    query_missing_GO_terms(created_entries)
    create_hierarchies(created_entries)

    notify_tasks.set_progress.apply_async((exp_id, max_progress, max_progress))

    return protein_map
