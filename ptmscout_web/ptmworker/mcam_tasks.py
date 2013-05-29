from ptmworker.helpers import upload_helpers
import celery
import zipfile
from ptmscout.config import settings
import os
from ptmworker import notify_tasks

def calculateEnrichment(experiment_id):
    pass

def enrichBool(output_filename, cluster_sets, enrichment, pvalue_cutoff):
    with open(output_filename, 'w') as ebf:
        for i, category in enumerate(sorted(enrichment.keys())):
            sorted_labels = sorted( enrichment[category]['labels'] )
            
            ebf.write("temp.type = '%s';\n" % (category))
            ebf.write("temp.labels = {'%s'};\n" % ("' '".join(sorted_labels)))
            
            bools = ";".join([ " ".join([ 1 if enrichment[category][label][cluster_set] < pvalue_cutoff else 0 for label in sorted_labels ]) for cluster_set in cluster_sets ])
            ebf.write("temp.bool = [%s];\n" % (bools))
            
            ebf.write("enrichBool{%d} = temp;\n" % (i+1))
        ebf.write("clear temp;")
                    
    
def numStruct(output_filename, cluster_sets, enrichment, pvalue_cutoff):
    with open(output_filename, 'w') as ebf:
        
        sorted_categories = sorted(enrichment.keys())
        ebf.write("numStruct.type = 'all';")
        ebf.write("numStruct.labels = {'%s'};" % ("' '".join( sorted_categories )))
        
        numStruct = ";".join([ " ".join([ sum([ 1 if enrichment[category][label][cluster_set] < pvalue_cutoff else 0 for label in sorted( enrichment[category]['labels'] ) ]) for category in sorted_categories ]) for cluster_set in cluster_sets ])
        
        ebf.write("numStruct.sum = [%s]" % (numStruct))

def zip(output_path, files):
    zf = zipfile.ZipFile(output_path, 'w')
    for f in files:
        zf.write(f)
    zf.close()

@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.logged_task
def run_mcam_analysis(output_filename, experiment_id, job_id):
    root_path = os.path.join(settings.ptmscout_path, settings.mcam_file_path, output_filename)
    os.mkdir(root_path)
    
    enrichment = calculateEnrichment(experiment_id)
    
    enrichBoolFilename = os.path.join(root_path, "enrichBool.m") 
    numStructFilename = os.path.join(root_path, "numStruct.m") 
    
    enrichBool(enrichBoolFilename, enrichment)
    numStruct(numStructFilename, enrichment)
    
    zippath = "%s.zip" % (root_path)
    zip(zippath, [enrichBoolFilename, numStruct])
    
    
    notify_tasks.finalize_mcam_export_job(job_id)