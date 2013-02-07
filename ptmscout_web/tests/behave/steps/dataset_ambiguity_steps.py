from behave import *
import bot

@given(u'a user chooses to assign default accessions for an experiment')
def assign_default_accessions(context):
    context.active_user.login()
    result = context.ptmscoutapp.get("/experiments/28/ambiguity?defaults=true")

    result = result.form.submit()
    result = result.follow()
    context.result = result.follow()

@when(u'the user submits the extension experiment for loading')
def submit_experiment(context):
    context.result = context.result.form.submit().follow()
    context.result = context.result.form.submit().follow()

    context.form = context.result.form
    context.form.set('terms_of_use', "yes")

@then(u'the user should be sent an email with a link to the extension experiment which contains')
@patch('transaction.abort')
@patch('transaction.commit')
@patch('ptmscout.utils.mail.celery_send_mail')
def experiment_uploaded_check_email(context, patch_mail, patch_commit, patch_abort):
    patch_commit.side_effect = bot.session_flush
    patch_abort.side_effect = bot.log_abort

    context.result = context.form.submit()

    exp_title = '[Default Assignments] Effects of HER2 overexpression on cell signaling networks governing proliferation and migration.'
    bot.check_experiment_loaded(context, exp_title, patch_mail)
