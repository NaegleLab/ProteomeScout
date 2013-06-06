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


def calculateEnrichment(scansite_cutoff, domain_cutoff, annotation_set_id, exp_id, user_id, job_id):
    measurements = modifications.getMeasuredPeptidesByExperiment(exp_id, secure=False)
    exp = experiment.getExperimentById(exp_id, secure=False)
    ptm_user = user.getUserById(user_id)
    
    annotation_types = dataset_explorer_view.compute_annotations(annotation_set_id, exp, ptm_user, measurements)
    
    cluster_sets = {}
    
    for ann_name in annotation_types:
        if annotation_types[ann_name] == 'cluster':
            cluster_sets[ann_name] = set()
    
    enrichment = {}
    
    for cluster_set in cluster_sets.keys():
        for ms in measurements:
            if cluster_set in ms.annotations:
                cluster_sets[cluster_set].add( ms.annotations[cluster_set] )
            
    for cluster_set in cluster_sets.keys():
        cluster_sets[cluster_set] = sorted( list( cluster_sets[cluster_set] ) ) 
        
    notify_tasks.set_job_stage.apply_async((job_id, 'Calculating Enrichment', len(cluster_sets)))
    
    test_cnt = defaultdict(lambda: 0)
    
    i = 0
    for cluster_set in cluster_sets.keys():
        for clabel in cluster_sets[cluster_set]:
            foreground = [ ms for ms in measurements if cluster_set in ms.annotations and ms.annotations[cluster_set] == clabel ]
            enrichment[(cluster_set, clabel)] = \
                dataset_explorer_view.calculate_feature_enrichment(foreground, measurements, annotation_types, 
                                                                   required_occurences=2, scansite_cutoff=scansite_cutoff, domain_cutoff=domain_cutoff )
            
            motif_cnt, motif_results = motif.calculate_motif_enrichment(foreground, measurements)
            enrichment[(cluster_set, clabel)] += [('motif', pep, pv) for pep, _fg, _bg, pv in motif_results]
            
            test_cnt['motif'] += motif_cnt
            
            
        i+=1
        notify_tasks.set_job_progress.apply_async((job_id, i, len(cluster_sets)))
        
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
            test_cnt[source] += 1
        
        nenrichment[(cluster_set, clabel)] = by_category
        
    return cluster_sets, nenrichment, enrichment_categories, test_cnt

@decorators.pushdir(os.path.join(settings.ptmscout_path, settings.mcam_file_path))
def enrichBool(output_filename, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff):
    
    with open(output_filename, 'w') as ebf:
        for i, category in enumerate(sorted(enrichment_categories.keys())):
            sorted_labels = enrichment_categories[category]
            
            ebf.write("temp.type = '%s';\n" % (category))
            ebf.write("temp.labels = {'%s'};\n" % ("' '".join(sorted_labels)))
            
            bools = ";".join([ 
                              " ".join([
                                        str( sum( 1 if label in enrichment[(cluster_set, clabel)][category] and enrichment[(cluster_set, clabel)][category][label] <= pvalue_cutoff else 0 for clabel in cluster_sets[cluster_set] ) )
                                            for label in sorted_labels ])
                                        for cluster_set in cluster_sets ])
            ebf.write("temp.bool = [%s];\n" % (bools))
            
            ebf.write("enrichBool{%d} = temp;\n" % (i+1))
        ebf.write("clear temp;\n")
                    
@decorators.pushdir(os.path.join(settings.ptmscout_path, settings.mcam_file_path))
def numStruct(output_filename, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff):
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
                                    for cluster_set in cluster_sets ])
        
        ebf.write("numStruct.sum = [%s];\n" % (numStruct))

@decorators.pushdir(os.path.join(settings.ptmscout_path, settings.mcam_file_path))
def archive(output_file, files):
    zf = zipfile.ZipFile(output_file, 'w')
    for f in files:
        zf.write(f)
    zf.close()

@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.logged_task
def run_mcam_analysis(output_filename, scansite_cutoff, domain_cutoff, pvalue_cutoff, correction_type, annotation_set_id, experiment_id, user_id, job_id):
    root_path = os.path.join(settings.ptmscout_path, settings.mcam_file_path, output_filename)
    os.mkdir(root_path)
    
    notify_tasks.set_job_status.apply_async((job_id, 'running'))
    
    cluster_sets, enrichment, enrichment_categories, test_cnt = calculateEnrichment(scansite_cutoff, domain_cutoff, annotation_set_id, experiment_id, user_id, job_id)
    
    if correction_type == 'fdr':
        fdrCorrection(pvalue_cutoff, enrichment, cluster_sets, enrichment_categories)
    if correction_type == 'bon':
        bonCorrection(enrichment, cluster_sets, enrichment_categories, test_cnt)
    
    enrichBoolFilename = os.path.join(output_filename, "loadEnrichBool.m") 
    numStructFilename = os.path.join(output_filename, "loadNumStruct.m") 
    
    notify_tasks.set_job_stage.apply_async((job_id, 'Writing Output', 0))
    
    enrichBool(enrichBoolFilename, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff)
    numStruct(numStructFilename, cluster_sets, enrichment, enrichment_categories, pvalue_cutoff)
    
    zipfilename = "%s.zip" % (output_filename)
    archive(zipfilename, [enrichBoolFilename, numStructFilename])
    
    
    notify_tasks.finalize_mcam_export_job.apply_async((job_id,))
