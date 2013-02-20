from ptmscout.config import settings, strings
from ptmscout.database import user
from ptmscout.database.user import NoSuchUser
from ptmscout.utils import webutils, mail, forms
from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

def create_schema(request):
    schema = forms.FormSchema()
    
    schema.add_text_field('username', "Username", 30)
    schema.add_password_field('pass1', "Password", maxlen=20)
    schema.add_password_field('pass2', "Repeat Password", maxlen=20)

    schema.add_text_field('name', "Name", 50)
    schema.add_text_field('email', "E-mail Address", 50)
    schema.add_text_field('institution', "Institution", 100)    

    schema.set_required_field('username')
    schema.set_required_field('pass1')
    schema.set_field_required_condition('pass2', 'pass1', forms.field_not_empty_test)
    schema.set_required_field('name')
    schema.set_required_field('email')
    schema.set_required_field('institution')
    
    schema.parse_fields(request)
    
    return schema

def validate_username(field_name, field_value, schema):
    try:
        _ptm_user = user.getUserByUsername(field_value)
        return strings.failure_reason_username_inuse
    except NoSuchUser:
        pass

def validate_password_length(field_name, field_value, schema):
    if len(field_value) < settings.MINIMUM_PASSWORD_LENGTH:
        return strings.failure_reason_password_too_short % settings.MINIMUM_PASSWORD_LENGTH
    
def validate_password_match(field_name, field_value, schema):
    pass1 = schema.get_form_value('pass1')
    if field_value != pass1:
        return strings.failure_reason_new_passwords_not_matching
        
def validate_email(field_name, field_value, schema):
    email_valid, domain_allowed = mail.email_is_valid(field_value)
    if not email_valid:
        return strings.failure_reason_email_not_valid
    if not domain_allowed:
        return strings.failure_reason_email_not_allowed
    

def create_validator(schema):
    validator = forms.FormValidator(schema)
    
    validator.add_alternate_validator('username', validate_username)
    validator.add_alternate_validator('pass1', validate_password_length)
    validator.add_alternate_validator('pass2', validate_password_match)
    validator.add_alternate_validator('email', validate_email)
    
    return validator


@view_config(route_name='register',renderer='ptmscout:templates/accounts/user_registration.pt')
def user_registration_view(request):
    submitted = webutils.post(request, 'submitted', "false") == "true"
    
    schema = create_schema(request)
    errors = []
    
    if submitted:
        errors = create_validator(schema).validate()
        
        if len(errors) == 0:
            __process_registration(request)
            return HTTPFound(request.application_url + "/registration_success")
    
    return {'errors': errors,
            'formrenderer': forms.FormRenderer(schema),
            'pageTitle': strings.user_registration_page_title}
    
@view_config(route_name='registration_success', renderer='ptmscout:templates/info/information.pt')
def user_registration_success(request):
    return {
        'pageTitle': strings.user_registration_page_title,
        'header': strings.user_registration_success_header,
        'message': strings.user_registration_success_message,
        'redirect': request.application_url + "/login"}


def __process_registration(request):
    username = webutils.post(request, 'username', "").strip()
    pass1    = webutils.post(request, 'pass1', "").strip()
    name      = webutils.post(request, 'name', "").strip()
    email     = webutils.post(request, 'email', "").strip()
    institute = webutils.post(request, 'institution', "").strip()
    
    ptm_user = user.User(username, name, email, institute)
    ptm_user.createUser(pass1)
    ptm_user.processInvitations()
    ptm_user.saveUser()
    
    mail.send_automail_message(request, 
    [email], 
    strings.user_registration_email_subject, 
    strings.user_registration_email_message % (ptm_user.name, request.application_url, ptm_user.username, ptm_user.activation_token))
