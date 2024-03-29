from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.utils import forms, webutils, decorators
from ptmscout.database import upload, experiment
from pyramid.httpexceptions import HTTPFound
from ptmscout.utils import wizard

MAX_VALUES=100
CONDITION_TYPES = \
    [('cell', "Cell Type"),('tissue',"Tissue Type"),
     ('drug', "Drug"), ('stimulus', "Stimulus"), ('environment', "Environmental Condition")]

def save_form_data(experiment, schema, added_fields):
    experiment.conditions = []
    
    for i in added_fields:
        t = schema.get_form_value('%d_type' % (i))
        v = schema.get_form_value('%d_value' % (i))
        experiment.addExperimentCondition(t, v)
        
    experiment.saveExperiment()
    
def parse_fields_from_experiment(exp, schema):
    for i, cond in enumerate(exp.conditions):
        schema.set_field('%d_type' % (i), cond.type)
        schema.set_field('%d_value' % (i), cond.value)

def get_form_schema(exp, parent_exp, request):
    submitted = webutils.post(request, 'submitted', False) == "true"
    schema = forms.FormSchema()
    
    for i in xrange(0, MAX_VALUES):
        type_field = '%d_type' % (i)
        value_field = '%d_value' % (i)
        schema.add_select_field(type_field, "Parameter Type %d" % (i+1), CONDITION_TYPES)
        schema.add_text_field(value_field, "Parameter Value %d" % (i+1), 50, 51)
        schema.set_field_required_condition(value_field, type_field, lambda type_value: type_value != '' and type_value != None)
    
    schema.parse_fields(request)
    if not submitted:
        if len(exp.conditions) == 0 and parent_exp != None:
            parse_fields_from_experiment(parent_exp, schema)
        else:
            parse_fields_from_experiment(exp, schema)
    
    added_fields = set()
    
    for i in xrange(0, MAX_VALUES):
        type_field = '%d_type' % (i)
        value_field = '%d_value' % (i)
        
        if schema.field_was_attempted(type_field) or schema.field_was_attempted(value_field):
            added_fields.add(i)
    
    return schema, added_fields

def mark_status(session):
    session.stage='confirm'
    session.save()


def create_nav_wizard(request, session):
    navigation = wizard.WizardNavigation(request)

    navigation.add_page('upload_config', "Configure Datafile", True, id=session.id)
    navigation.add_page('upload_metadata', "Add Metadata", True, id=session.id)
    navigation.add_page('upload_conditions', "Describe Conditions", True, id=session.id)
    navigation.add_page('upload_confirm', "Confirm Upload", False, id=session.id)
    navigation.set_page('upload_conditions')

    return navigation

def upload_conditions(request, session):
    submitted = webutils.post(request, 'submitted', False) == "true"
    
    exp = experiment.getExperimentById(session.experiment_id, request.user, False)
    parent_exp = None

    if session.parent_experiment != None:
        parent_exp = experiment.getExperimentById(session.parent_experiment, request.user, False)
    
    navigation = create_nav_wizard(request, session)
    schema, added_fields = get_form_schema(exp, parent_exp, request)
    renderer = forms.FormRenderer(schema)
    
    errors = []
    
    if submitted:
        errors = forms.FormValidator(schema).validate()
        if len(errors) == 0:
            save_form_data(exp, schema, added_fields)
            mark_status(session)
            return HTTPFound(request.application_url + "/upload/%d/confirm" % (session.id))
    
    return {'formrenderer':renderer,
            'navigation':navigation,
            'pageTitle': strings.experiment_upload_conditions_page_title,
            'header': strings.experiment_upload_conditions_page_title,
            'errors':errors,
            'session_id': session.id,
            'added_fields':added_fields,
            'MAX_FIELDS':MAX_VALUES
            }

@view_config(route_name='upload_conditions', renderer='ptmscout:/templates/upload/upload_conditions.pt')
@decorators.get_session('id','experiment')
def upload_conditions_view(context, request, session):
    return upload_conditions(request, session)
