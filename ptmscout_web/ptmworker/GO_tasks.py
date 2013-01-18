import celery
import logging
from ptmworker import notify_tasks
from ptmworker.helpers import upload_helpers, quickgo_tools
from ptmscout.database import protein

log = logging.getLogger('ptmscout')

MAX_QUICKGO_BATCH_SIZE = 500

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


def create_missing_GO_annotations(GO_annotation_map, accession_map, protein_map):
    created_go_entries = []
    created, assigned = 0, 0
    
    complete_go_terms = set()
    for acc in GO_annotation_map:
        for goId in GO_annotation_map[acc]:
            complete_go_terms.add(goId)

    missing_terms = set()
    unique_terms = set()
    for acc in accession_map:
        go_annotations = GO_annotation_map[acc]

        for protein_id in accession_map[acc]:
            prot = protein_map[protein_id]

            for goId in go_annotations:
                unique_terms.add(goId)

                if not prot.hasGoTerm(goId):
                    dateAdded = go_annotations[goId]
                    go_term, entry = upload_helpers.get_go_annotation(goId, complete_go_terms, missing_terms)

                    prot.addGoTerm(go_term, dateAdded)

                    if entry:
                        created_go_entries.append(entry)
                        created+=1

                    assigned+=1

    for pid in protein_map:
        protein_map[pid].saveNoFlush()

    missing_terms = list(missing_terms)
    upload_helpers.query_missing_GO_terms(missing_terms, complete_go_terms)

    log.info("Returned %d terms, assigned %d terms, created %d terms", len(unique_terms), assigned, created)
    return created_go_entries


def build_accession_map(protein_ids):
    valid_types = set(['refseq','uniprot','swissprot','gene_synonym'])

    protein_accession_map = {}
    protein_map = {}
    for protein_id in protein_ids:
        prot = protein.getProteinById(protein_id)
        protein_map[protein_id] = prot

        for accession in prot.accessions:
            if accession.type in valid_types:
                pids = protein_accession_map.get(accession.value, set())
                pids.add(protein_id)
                protein_accession_map[accession.value] = pids

    return protein_accession_map, protein_map


@celery.task
@upload_helpers.transaction_task
def import_go_terms(protein_result, exp_id):
    protein_map, new_protein_ids = protein_result

    protein_accession_multimap, protein_id_map = build_accession_map(new_protein_ids.values())
    GO_map = {}
    GO_annotation_task_args = upload_helpers.create_chunked_tasks(protein_accession_multimap, MAX_QUICKGO_BATCH_SIZE)
    max_progress = len(GO_annotation_task_args) + 1

    upload_helpers.store_stage_input(exp_id, 'GO terms', protein_result)
    notify_tasks.set_loading_stage.apply_async((exp_id, 'GO terms', max_progress))

    i=0
    for protein_accessions in GO_annotation_task_args:
        go_terms, _ = quickgo_tools.batch_get_GO_annotations(protein_accessions)
        GO_map.update( go_terms )
        i+=1
        notify_tasks.set_progress.apply_async((exp_id, i, max_progress))

    created_entries = create_missing_GO_annotations(GO_map, protein_accession_multimap, protein_id_map)
    create_hierarchy_for_missing_GO_annotations(created_entries)

    notify_tasks.set_progress.apply_async((exp_id, max_progress, max_progress))

    return protein_map
