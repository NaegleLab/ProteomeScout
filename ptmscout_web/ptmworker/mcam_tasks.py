from ptmworker.helpers import upload_helpers
import celery
import zipfile
from ptmscout.config import settings
import os
from ptmworker import notify_tasks
from ptmscout.database import experiment, modifications, user
from ptmscout.views.dataset import dataset_explorer_view
from collections import defaultdict
from ptmscout.utils import decorators, motif
import cProfile

MCAM_SUBPROCESSES = 5

def fdrCorrection(pvalue_cutoff, enrichment, cluster_sets, enrichment_categories):
    test_category_pvalues = {}

    for category in enrichment_categories:
        test_category_pvalues[category] = []
        for label in enrichment_categories[category]:
            for cluster_set in cluster_sets:
                for clabel in cluster_sets[cluster_set]:
                    if label in enrichment[(cluster_set, clabel)][category]:
                        test_category_pvalues[category].append( enrichment[(cluster_set, clabel)][category][label] )
                    
    test_category_cutoffs = {}
    for category in list(test_category_pvalues.keys()):
        sorted_pvalues = sorted( test_category_pvalues[category] )
        
        m = len(sorted_pvalues)
        bh_i = m - 1
        for i, p in enumerate(sorted_pvalues):
            if p < (i / float(m)) * pvalue_cutoff:
                bh_i = i
                break
        
        test_category_cutoffs[category] = sorted_pvalues[bh_i]

    for cluster_set in cluster_sets:
        for clabel in cluster_sets[cluster_set]:
            for category in enrichment_categories:
                for label in enrichment_categories[category]:
                    if label in enrichment[(cluster_set, clabel)][category] and \
                            enrichment[(cluster_set, clabel)][category][label] > test_category_cutoffs[category]:
                        enrichment[(cluster_set, clabel)][category][label] = 1.0

def bonCorrection(enrichment, cluster_sets, enrichment_categories, test_cnt):
    for cluster_set in cluster_sets:
        for clabel in cluster_sets[cluster_set]:
            for category in enrichment_categories:
                for label in enrichment_categories[category]:
                    if label in enrichment[(cluster_set, clabel)][category]: 
                        enrichment[(cluster_set, clabel)][category][label] *= test_cnt[category]

@celery.task
@upload_helpers.transaction_task
def calculate_feature_enrichment(cluster_set_list, cluster_sets, annotation_types, scansite_cutoff, domain_cutoff, exp_id, job_id):
    measurements = modifications.getMeasuredPeptidesByExperiment(exp_id, secure=False)
    test_cnt = defaultdict(lambda: 0)
    label_cache = {'proteins': defaultdict(dict), 'peptides': defaultdict(dict)}
    enrichment = {}

    for cluster_set in cluster_set_list:
        for clabel in cluster_sets[cluster_set]:
            foreground = [ ms for ms in measurements if cluster_set in ms.annotations and ms.annotations[cluster_set] == clabel ]
            enrichment[(cluster_set, clabel)] = \
                dataset_explorer_view.calculate_feature_enrichment(foreground, measurements, annotation_types,
                                                                   required_occurences=2,
                                                                   scansite_cutoff=scansite_cutoff,
                                                                   domain_cutoff=domain_cutoff,
                                                                   cache_table = label_cache)
            
            motif_cnt, motif_results = motif.calculate_motif_enrichment(foreground, measurements)
            enrichment[(cluster_set, clabel)] += [('motif', pep, pv) for pep, _fg, _bg, pv in motif_results]
            
            test_cnt['motif'] += motif_cnt
            notify_tasks.increment_job_progress.apply_async((job_id,))

    return enrichment, test_cnt

def aggregate_results(result):
    aggregate_enrichment = {}
    aggregate_test_cnt = defaultdict(lambda: 0)

    for enrichment, test_cnt in result.join():
        aggregate_enrichment.update(enrichment)
        for category in test_cnt:
            aggregate_test_cnt[category] += test_cnt[category]

    return aggregate_enrichment, aggregate_test_cnt

def calculateEnrichment(scansite_cutoff, domain_cutoff, annotation_set_id, exp_id, user_id, job_id):
    measurements = modifications.getMeasuredPeptidesByExperiment(exp_id, secure=False)
    exp = experiment.getExperimentById(exp_id, secure=False)
    ptm_user = user.getUserById(user_id)
    
    annotation_types, annotation_order = dataset_explorer_view.compute_annotations(annotation_set_id, exp, ptm_user, measurements)
    
    cluster_sets = {}
    
    for ann_name in annotation_types:
        if annotation_types[ann_name] == 'cluster':
            cluster_sets[ann_name] = set()
    
    sorted_clusters = sorted( cluster_sets.keys(), key=lambda ann_name: annotation_order[ann_name] )
    
    for cluster_set in cluster_sets.keys():
        for ms in measurements:
            if cluster_set in ms.annotations:
                cluster_sets[cluster_set].add( ms.annotations[cluster_set] )
            
    max_progress = 0
    for cluster_set in cluster_sets.keys():
        cluster_sets[cluster_set] = sorted( list( cluster_sets[cluster_set] ) )
        max_progress += len(cluster_sets[cluster_set])

    notify_tasks.set_job_stage.apply_async((job_id, 'Calculating Enrichment', max_progress))
    
    cluster_chunks = upload_helpers.create_chunked_tasks( cluster_sets.keys(), len(cluster_sets.keys()) / MCAM_SUBPROCESSES )
    
    cluster_enrichment_jobs = celery.group([ calculate_feature_enrichment.s(cluster_chunk, cluster_sets, annotation_types, scansite_cutoff, domain_cutoff, exp_id, job_id) for cluster_chunk in cluster_chunks ])

    result = cluster_enrichment_jobs.apply_async()
    enrichment, test_cnt = aggregate_results(result)

    enrichment_categories = defaultdict(set)
    for (cluster_set, clabel) in enrichment:
        for source, label, p_value in enrichment[(cluster_set, clabel)]:
            enrichment_categories[source].add(label)
            
    for source in enrichment_categories.keys():
        enrichment_categories[source] = sorted( list( enrichment_categories[source] ) )
        
    nenrichment = {}
    for (cluster_set, clabel) in enrichment.keys():
        by_category = defaultdict(dict)

        table = enrichment[(cluster_set, clabel)]
        for source, label, p_value in table:
            by_category[source][label] = p_value
            if source != 'motif':
                test_cnt[source] += 1
        
        nenrichment[(cluster_set, clabel)] = by_category
        
    return sorted_clusters, cluster_sets, nenrichment, enrichment_categories, test_cnt

@decorators.pushdir(os.path.join(settings.ptmscout_path, settings.mcam_file_path))
def enrichBool(output_filename, sorted_cluster_sets, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff):
    
    with open(output_filename, 'w') as ebf:
        for i, category in enumerate(sorted(enrichment_categories.keys())):
            sorted_labels = enrichment_categories[category]
            
            ebf.write("temp.type = '%s';\n" % (category))
            ebf.write("temp.labels = {'%s'};\n" % ("' '".join(sorted_labels)))
            
            bools = ";".join([ 
                              " ".join([
                                        str( sum( 1 if label in enrichment[(cluster_set, clabel)][category] and enrichment[(cluster_set, clabel)][category][label] <= pvalue_cutoff else 0 for clabel in cluster_sets[cluster_set] ) )
                                            for label in sorted_labels ])
                                        for cluster_set in sorted_cluster_sets ])
            ebf.write("temp.bool = [%s];\n" % (bools))
            
            ebf.write("enrichBool{%d} = temp;\n" % (i+1))
        ebf.write("clear temp;\n")
                    
@decorators.pushdir(os.path.join(settings.ptmscout_path, settings.mcam_file_path))
def numStruct(output_filename, sorted_cluster_sets, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff):
    with open(output_filename, 'w') as ebf:
        
        sorted_categories = sorted(enrichment_categories.keys())
        ebf.write("numStruct.type = 'all';\n")
        ebf.write("numStruct.labels = {'%s'};\n" % ("' '".join( sorted_categories )))
        
        numStruct = ";".join([
                              " ".join([
                                        str(
                                            sum([
                                                max([ 1 if label in enrichment[(cluster_set, clabel)][category] and enrichment[(cluster_set,clabel)][category][label] <= pvalue_cutoff else 0
                                                        for clabel in cluster_sets[cluster_set] ])
                                                 for label in enrichment_categories[category] ]) )
                                        for category in sorted_categories ])
                                    for cluster_set in sorted_cluster_sets ])
        
        ebf.write("numStruct.sum = [%s];\n" % (numStruct))

@decorators.pushdir(os.path.join(settings.ptmscout_path, settings.mcam_file_path))
def archive(output_file, files):
    zf = zipfile.ZipFile(output_file, 'w')
    for f in files:
        zf.write(f)
    zf.close()

def filterUnenriched(cluster_sets, enrichment, enrichment_categories, pvalue_cutoff):

    for category in enrichment_categories:
        remove_labels = set()
        for label in enrichment_categories[category]:
            has_enriched = False
            for cluster_set in cluster_sets:
                for clabel in cluster_sets[cluster_set]:
                    has_enriched |= ( label in enrichment[(cluster_set, clabel)][category] and enrichment[(cluster_set,clabel)][category][label] <= pvalue_cutoff )
            if not has_enriched:
                remove_labels.add(label)

        for label in remove_labels:
            enrichment_categories[category].remove(label)

@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.transaction_task
def run_mcam_analysis(output_filename, scansite_cutoff, domain_cutoff, pvalue_cutoff, correction_type, annotation_set_id, experiment_id, user_id, job_id):
    root_path = os.path.join(settings.ptmscout_path, settings.mcam_file_path, output_filename)
    os.mkdir(root_path)
    
    notify_tasks.set_job_status.apply_async((job_id, 'started'))
    
    sorted_cluster_sets, cluster_sets, enrichment, enrichment_categories, test_cnt = calculateEnrichment(scansite_cutoff, domain_cutoff, annotation_set_id, experiment_id, user_id, job_id)
    
    if correction_type == 'fdr':
        fdrCorrection(pvalue_cutoff, enrichment, cluster_sets, enrichment_categories)
    if correction_type == 'bon':
        bonCorrection(enrichment, cluster_sets, enrichment_categories, test_cnt)
   
    filterUnenriched(cluster_sets, enrichment, enrichment_categories, pvalue_cutoff)

    enrichBoolFilename = os.path.join(output_filename, "loadEnrichBool.m") 
    numStructFilename = os.path.join(output_filename, "loadNumStruct.m") 
    
    notify_tasks.set_job_stage.apply_async((job_id, 'Writing Output', 0))
    
    enrichBool(enrichBoolFilename, sorted_cluster_sets, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff)
    numStruct(numStructFilename, sorted_cluster_sets, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff)
    
    zipfilename = "%s.zip" % (output_filename)
    archive(zipfilename, [enrichBoolFilename, numStructFilename])
    
    
    notify_tasks.finalize_mcam_export_job.apply_async((job_id,))
