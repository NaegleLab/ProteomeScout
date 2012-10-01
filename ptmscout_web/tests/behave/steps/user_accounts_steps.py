from bot import Bot
from mock import patch
import re
from ptmscout.config import strings

@given(u'I want to create an account and I have a .edu email address')
def setup_variables(context):
    context.active_user = Bot(context.ptmscoutapp)

@given(u'I have registered for an account')
def register_account(context):
    context.active_user = Bot(context.ptmscoutapp)
    context.active_user.register()
    

@given(u'I have an existing account')
def register_and_activate(context):
    context.active_user = Bot(context.ptmscoutapp)
    context.active_user.register_and_login()
    context.active_user.logout()
    
    
@when(u'I enter a my user information')
def fill_out_registration_form(context):
    context.result = context.ptmscoutapp.get('/register', status=200)
    context.form = context.result.form
    
    context.form.set('username', context.active_user.username)
    context.form.set('pass1', context.active_user.password)
    context.form.set('pass2', context.active_user.password)
    context.form.set('name', context.active_user.fullname)
    context.form.set('email', context.active_user.email)
    context.form.set('institution', context.active_user.institution)

@when(u'I follow the link in my e-mail to activate my account')
def follow_activation_link(context):
    response = context.ptmscoutapp.get(context.active_user.activation_link, status=200)
    response.mustcontain(strings.account_activation_success_header)

@when(u'I press "I forgot my password" on the login page and enter my email address')
def choose_forgot_password(context):
    context.result = context.ptmscoutapp.get('/login', status=200)
    context.result = context.result.click(href='http://localhost/forgot_password')
    
    context.form = context.result.form
    context.form.set("email", context.active_user.email)
    
    
@then(u'I should be sent an e-mail on how to activate my account')
@patch('ptmscout.utils.mail.send_automail_message')
def intercept_mail(context, patch_mail):
    context.result = context.form.submit()
    patch_mail.assert_called_once()

@then(u'I should be able to log in to my account')
def login(context):
    context.result = context.ptmscoutapp.get('/login', status=200)
    context.form = context.result.form
    
    context.form.set("username", context.active_user.username)
    context.form.set("password", context.active_user.password)
    context.result = context.form.submit()
    
    context.result.mustcontain(strings.login_page_success_message)

@then(u'I should be sent a temporary password to login')
@patch('ptmscout.utils.mail.send_automail_message')
def reset_password_received(context, patch_mail):
    context.result = context.form.submit()
    patch_mail.assert_called_once()
    
    regex = re.compile("Username: ([a-zA-Z0-9\_\-]+)\\\\nPassword: ([a-z0-9A-Z]+)\\\\n")
    regex_result = regex.search(str(patch_mail.call_args))
    
    new_username = regex_result.group(1)
    new_password = regex_result.group(2)
    
    context.active_user.username = new_username
    context.active_user.password = new_password
    context.active_user.login()
    