from pyramid.view import view_config
from ptmscout.config import strings, settings
from ptmscout.utils import forms, decorators
from pyramid.httpexceptions import HTTPFound, HTTPNotFound
from ptmscout.database import annotations, jobs
from ptmworker import mcam_tasks
import time
from pyramid.response import FileResponse
import os

@view_config(route_name='mcam_download', permission='private')
@decorators.get_experiment('id')
def download_mcam_file(context, request, experiment):
    fname = "%s.mcam.%d.%s.zip" % (request.matchdict['id'], request.user.id, request.matchdict['mcam_id'])
    
    fpath = os.path.join(settings.ptmscout_path, settings.mcam_file_path, fname)
    if not os.path.exists(fpath):
        raise HTTPNotFound()

    response = FileResponse(fpath, request)
    response.content_type = 'application/zip'
    response.content_disposition = 'attachment; filename="%s"' % (fname)
    return response

@view_config(route_name='mcam_confirm', renderer='ptmscout:/templates/info/information.pt', permission='private')
def upload_already_started_view(request):
    return {'pageTitle': strings.mcam_enrichment_started_page_title,
            'header': strings.mcam_enrichment_started_page_title,
            'message': strings.mcam_enrichment_started_message % (request.route_url('my_experiments'))}

def get_annotation_sets(experiment, user):
    annotation_sets = []
    for annotation_set in annotations.getUserAnnotationSets(experiment.id, user):
        annotation_sets.append( (str(annotation_set.id), annotation_set.name) )
    return annotation_sets

def start_mcam_job(request, schema, experiment, user):
    mcam_id = int(time.time())
    mcam_filename_base = "%d.mcam.%d.%d" % (experiment.id, user.id, mcam_id)
    
    mcam_job = jobs.Job()
    mcam_job.name = "Run MCAM Enrichment for Experiment %d" % (experiment.id)
    mcam_job.user = user
    mcam_job.result_url = request.route_url('mcam_download', id=experiment.id, mcam_id=mcam_id)
    mcam_job.status_url = request.route_url('my_experiments')
    mcam_job.status = 'in queue'
    mcam_job.type = 'mcam_enrichment'
    
    mcam_job.save()
    
    annotation_set_id = int( schema.get_form_value('annotationset') ) 
    scansite_cutoff = float( schema.get_form_value('scansitecutoff') ) 
    domain_cutoff = float( schema.get_form_value('domaincutoff') )
    pvalue_cutoff = float( schema.get_form_value('alpha') )
    correction_type = schema.get_form_value('correction')
    
    mcam_tasks.run_mcam_analysis.apply_async((mcam_filename_base, scansite_cutoff, domain_cutoff, pvalue_cutoff, correction_type, annotation_set_id, experiment.id, user.id, mcam_job.id))

def create_schema(request, annotation_sets):
    schema = forms.FormSchema()
    
    schema.add_select_field('annotationset', "Annotation Set", annotation_sets)
    schema.add_radio_field('correction', "Correction Type", [('fdr',"BH FDR"),('bon',"Bonferroni")], default='fdr')

    schema.add_decimal_field('scansitecutoff', "Scansite Cutoff", maxlen=55, default=str(settings.default_mcam_scansite_cutoff))
    schema.add_decimal_field('domaincutoff', "Domain Cutoff", maxlen=55, default=str(settings.default_mcam_domain_cutoff))
    schema.add_decimal_field('alpha', "Alpha P-value", maxlen=55, default=str(settings.default_mcam_p_cutoff))
    
    schema.set_required_field('annotationset')
    schema.set_required_field('correction')
    schema.set_required_field('domaincutoff')
    schema.set_required_field('scansitecutoff')
    schema.set_required_field('alpha')
    
    schema.parse_fields(request)
    
    return schema

@view_config(route_name='mcam_enrichment', request_method='POST', renderer='ptmscout:/templates/experiments/experiment_mcam.pt', permission='private', )
@decorators.get_experiment('id')
def view_mcam_enrichment_POST(context, request, experiment):
    annotation_sets = get_annotation_sets(experiment, request.user)
    schema = create_schema(request, annotation_sets)
    errors = forms.FormValidator(schema).validate()
    
    if len(errors) > 0:
        return {'pageTitle': strings.mcam_enrichment_page_title,
                'formrenderer': forms.FormRenderer(schema),
                'experiment':experiment,
                'errors':errors}
        
    start_mcam_job(request, schema, experiment, request.user)
    return HTTPFound(request.route_url('mcam_confirm', id=experiment.id))

@view_config(route_name='mcam_enrichment', request_method='GET', renderer='ptmscout:/templates/experiments/experiment_mcam.pt', permission='private', )
@decorators.get_experiment('id',types=['dataset','experiment'])
def view_mcam_enrichment_GET(context, request, experiment):
    annotation_sets = get_annotation_sets(experiment, request.user)
    schema = create_schema(request, annotation_sets)
    
    return {'pageTitle': strings.mcam_enrichment_page_title,
            'formrenderer': forms.FormRenderer(schema),
            'experiment':experiment,
            'errors':[]}
