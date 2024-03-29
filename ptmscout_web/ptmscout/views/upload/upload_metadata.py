from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.utils import webutils, forms, decorators, wizard
from pyramid.httpexceptions import HTTPFound
from ptmscout.database import experiment, upload

def write_experiment_properties(nexp, session, schema, current_user):
    nexp.name = schema.get_form_value('experiment_name')
    nexp.description = schema.get_form_value('description')
    
    nexp.dataset = session.data_file
    
    nexp.published = 1 if schema.get_form_value('published') == "yes" else 0
    nexp.ambiguity = 1 if schema.get_form_value('ambiguous') == "yes" else 0
    nexp.URL = schema.get_form_value('URL')
    nexp.type = 'experiment'

    nexp.contact = None
    nexp.author = None
    nexp.journal = None
    nexp.volume = None
    nexp.page_end = None
    nexp.page_start = None
    nexp.publication_month = None
    nexp.publication_year = None
    nexp.errorLog = None
    nexp.PMID = None
    
    nexp.submitter_id = current_user.id
    
    if(session.load_type == 'extension'):
        nexp.experiment_id = session.parent_experiment
        
    if(schema.get_form_value('published') == "yes"):
        nexp.contact = schema.get_form_value('author_contact')
        nexp.author = schema.get_form_value('authors')
        
        nexp.volume = int(schema.get_form_value('volume'))
        nexp.page_start = schema.get_form_value('page_start')
        nexp.page_end = schema.get_form_value('page_end')
        
        nexp.publication_month = schema.get_form_value('publication_month')
        nexp.publication_year = int(schema.get_form_value('publication_year'))
        
        nexp.journal = schema.get_form_value('journal')
        if(schema.get_form_value('pmid') != ""):
            nexp.PMID = int(schema.get_form_value('pmid'))
    
    
    nexp.grantPermission(current_user, 'owner')
    nexp.saveExperiment()


def get_experiment_ref(session, current_user):
    if(session.experiment_id != None):
        nexp = experiment.getExperimentById(session.experiment_id, current_user, False)
    else: 
        nexp = experiment.Experiment()
        nexp.public = 0
        nexp.experiment_id = None
    
    return nexp
    
def mark_status(nexp, session):
    session.experiment_id = nexp.id
    session.stage = 'condition'
    session.save()

def create_schema(request):
    schema = forms.FormSchema()
    
    schema.add_text_field('experiment_name', "Experiment Name", width=55)
    schema.add_text_field('URL', "URL", width=55)
    schema.add_textarea_field('description', "Description", 42, 5)
    
    schema.add_select_field('published', 'Published', [('no',"No"), ('yes',"Yes")])
    schema.add_integer_field('pmid', 'PubMed Id', maxlen=10)
    
    
    schema.add_text_field('authors', 'Authors', width=55)
    schema.add_text_field('author_contact', 'Author Contact E-mail', width=55)
    months = ["january","february", "march", "april",
              "may","june","july","august","september",
              "october","november","december"]
    schema.add_select_field('publication_month', 'Publication Month', [(m, m.capitalize()) for m in months])
    schema.add_integer_field('publication_year', 'Publication Year', maxlen=4)
    schema.add_text_field('journal', 'Journal', width=55)
    schema.add_integer_field('volume', 'Volume', 10)
    schema.add_text_field('page_start', 'Page Start', maxlen=10, width=11)
    schema.add_text_field('page_end', 'Page End', maxlen=10, width=11)
    
    schema.add_select_field('ambiguous', 'Ambiguous', [('no',"No"), ('yes',"Yes")])
    
    schema.set_required_field('description')
    schema.set_required_field('experiment_name')
    schema.set_required_field('published')
    schema.set_required_field('ambiguous')
    
    schema.set_field_required_condition('author_contact', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('authors', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('journal', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('publication_year', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('publication_month', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('volume', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('page_start', 'published', forms.field_equals_test('yes'))
    schema.set_field_required_condition('page_end', 'published', forms.field_equals_test('yes'))
    
    schema.parse_fields(request)
    
    return schema
    
def populate_schema_from_experiment(schema, session, user):
    field_dict = None
    if session.experiment_id != None:
        exp = experiment.getExperimentById(session.experiment_id, user, False)
        field_dict = {
                     'pmid':exp.PMID,
                     'publication_year': exp.publication_year,
                     'publication_month': exp.publication_month,
                     'published': 'yes' if exp.published == 1 else 'no',
                     'ambiguous': 'yes' if exp.ambiguity == 1 else 'no',
                     'author_contact' : exp.contact,
                     'page_start': exp.page_start,
                     'page_end': exp.page_end,
                     'authors': exp.author,
                     'journal': exp.journal,
                     'volume': exp.volume,
                     'description': exp.description,
                     'experiment_name':exp.name,
                     'URL': exp.URL
                  }
    elif session.parent_experiment != None:
        exp = experiment.getExperimentById(session.parent_experiment, user, False)
        name = exp.name
        description = exp.description

        if 'extension' == session.load_type:
            name = "[%s] %s" % (session.change_name, exp.name)
            description = "%s -- %s" % (exp.description, session.change_description)

        field_dict = {
                     'pmid':exp.PMID,
                     'publication_year': exp.publication_year,
                     'publication_month': exp.publication_month,
                     'published': 'yes' if exp.published == 1 else 'no',
                     'ambiguous': 'yes' if exp.ambiguity == 1 else 'no',
                     'author_contact' : exp.contact,
                     'page_start': exp.page_start,
                     'page_end': exp.page_end,
                     'authors': exp.author,
                     'journal': exp.journal,
                     'volume': exp.volume,
                     'description': description,
                     'experiment_name': name,
                     'URL': exp.URL
                  }

        
    if field_dict != None:
        schema.parse_fields_dict(field_dict)

def create_nav_wizard(request, session):
    navigation = wizard.WizardNavigation(request)

    navigation.add_page('upload_config', "Configure Datafile", True, id=session.id)
    navigation.add_page('upload_metadata', "Add Metadata", True, id=session.id)
    navigation.add_page('upload_conditions', "Describe Conditions", False, id=session.id)
    navigation.add_page('upload_confirm', "Confirm Upload", False, id=session.id)
    navigation.set_page('upload_metadata')

    return navigation

def upload_metadata(request, session):
    submitted = webutils.post(request,'submitted',"false")
    
    schema = create_schema(request)
    errors = []
    
    if(submitted == "true"):
        errors = forms.FormValidator(schema).validate()
        
        if len(errors) == 0:        
            # need to get exp_file from upload session records
            nexp = get_experiment_ref(session, request.user)
            write_experiment_properties(nexp, session, schema, request.user)
            mark_status(nexp, session)
            return HTTPFound(request.application_url + "/upload/%d/conditions" % session.id)
    else:
        populate_schema_from_experiment(schema, session, request.user)
        
    return {'pageTitle': strings.upload_page_title,
            'navigation': create_nav_wizard(request, session),
            'errors':errors,
            'session_id': session.id,
            'formrenderer':forms.FormRenderer(schema)}

@view_config(route_name='upload_metadata', renderer='ptmscout:templates/upload/upload_metadata.pt', permission='private')
@decorators.get_session('id', 'experiment')
def upload_metadata_view(context, request, session):
    return upload_metadata(request, session)
