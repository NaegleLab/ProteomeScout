from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.utils import webutils, forms
from pyramid.httpexceptions import HTTPFound
from ptmscout.database import experiment, upload

def write_experiment_properties(nexp, session, schema, current_user):
    nexp.name = schema.get_form_value('experiment_name')
    nexp.description = schema.get_form_value('description')
    
    nexp.dataset = session.data_file
    nexp.status='preload'
    
    nexp.published = 1 if schema.get_form_value('published') == "yes" else 0
    nexp.ambiguity = 1 if schema.get_form_value('ambiguous') == "yes" else 0
    nexp.URL = schema.get_form_value('URL')
    nexp.export = 1

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
        nexp = experiment.getExperimentById(session.experiment_id, current_user)
    elif(session.load_type == 'reload' or session.load_type == 'append'):
        nexp = experiment.getExperimentById(session.parent_experiment, current_user)
    else: 
        nexp = experiment.Experiment()
        nexp.public = 0
        nexp.experiment_id = None
    
    return nexp
    
def mark_status(nexp, session):
    session.experiment_id = nexp.id
    session.stage = 'confirm'
    session.save()

def create_schema(request):
    schema = forms.FormSchema()
    
    schema.add_text_field('experiment_name', "Experiment Name", width=55)
    schema.add_text_field('URL', "URL", width=55)
    schema.add_textarea_field('description', "Description", 42, 5)
    
    schema.add_select_field('published', 'Published', [('no',"No"), ('yes',"Yes")])
    schema.add_numeric_field('pmid', 'PubMed Id', maxlen=10)
    
    
    schema.add_text_field('authors', 'Authors', width=55)
    schema.add_text_field('author_contact', 'Author Contact E-mail', width=55)
    months = ["january","february", "march", "april",
              "may","june","july","august","september",
              "october","november","december"]
    schema.add_select_field('publication_month', 'Publication Month', [(m, m.capitalize()) for m in months])
    schema.add_numeric_field('publication_year', 'Publication Year', maxlen=4)
    schema.add_text_field('journal', 'Journal', width=55)
    schema.add_numeric_field('volume', 'Volume', 10)
    schema.add_text_field('page_start', 'Page Start', maxlen=10, width=11)
    schema.add_text_field('page_end', 'Page End', maxlen=10, width=11)
    
    schema.add_select_field('ambiguous', 'Ambiguous', [('no',"No"), ('yes',"Yes")])
    schema.add_textarea_field('notes', "Notes", 42, 5)
    
    
    schema.set_required_field('description')
    schema.set_required_field('experiment_name')
    schema.set_required_field('published')
    schema.set_required_field('ambiguous')
    
    schema.set_field_required_condition('author_contact', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('authors', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('journal', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('publication_year', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('publication_month', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('volume', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('page_start', 'published', forms.field_not_empty_test)
    schema.set_field_required_condition('page_end', 'published', forms.field_not_empty_test)
    
    schema.parse_fields(request)
    
    return schema
    
def populate_schema_from_experiment(schema, session, user):
    exp = None
    if session.experiment_id != None:
        exp = experiment.getExperimentById(session.experiment_id, user)
    elif session.parent_experiment != None:
        exp = experiment.getExperimentById(session.parent_experiment, user)
        
    if exp != None:
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
                     'experiment_name':"",
                     'URL': exp.URL,
                     'notes':""
                  }
        schema.parse_fields_dict(field_dict)

@view_config(route_name='upload_metadata', renderer='ptmscout:templates/upload/upload_metadata.pt', permission='private')
def upload_metadata(request):
    submitted = webutils.post(request,'submitted',"false")
    session_id = int(request.matchdict['id'])
    session = upload.getSessionById(session_id, request.user)
    
    schema = create_schema(request)
    errors = []
    
    if(submitted == "true"):
        errors = forms.FormValidator(schema).validate()
        
        if len(errors) == 0:        
            # need to get exp_file from upload session records
            nexp = get_experiment_ref(session, request.user)
            write_experiment_properties(nexp, session, schema, request.user)
            mark_status(nexp, session)
            return HTTPFound(request.application_url + "/upload/%d/conditions" % session_id)
    else:
        populate_schema_from_experiment(schema, session, request.user)
        
    return {'pageTitle': strings.upload_page_title,
            'errors':errors,
            'formrenderer':forms.FormRenderer(schema)}
