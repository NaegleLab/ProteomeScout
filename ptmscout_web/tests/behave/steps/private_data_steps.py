from behave import *
from bot import Bot
from assertions import assertDoesNotContain, assertContains
from ptmscout.config import strings
from mock import patch
from ptmscout.database import permissions
import assertions
import helpers

@given(u'I have loaded a dataset and marked it private')
def create_owning_user_with_datasets(context):
    owner_user = Bot(context.ptmscoutapp)
    owner_user.register()
    owner_user.activate()
    exp_ids=[26,28]
    owner_user.acquire_experiments(exp_ids)
    
    context.owner_user = owner_user

@given(u'a user has been invited to view a private dataset')
def create_and_invite_unregistered_user(context):
    create_owning_user_with_datasets(context)
    invite_unregistered_user(context)
    check_invite_email_sent(context)

@when(u'I enter an email address for an unregistered user in "Share dataset"')
def invite_unregistered_user(context):
    context.owner_user.login()
    context.result = context.ptmscoutapp.get('/account/experiments/26/share', status=200)
    form = context.result.forms[0]
    
    context.invited_email = "invited@institute.edu"
    
    form.set('email', context.invited_email)
    context.result = form.submit()
    
    context.result = context.result.follow()
    
    context.form = context.result.forms[0]
    context.result.mustcontain(strings.user_invite_confirm % (context.invited_email))


@when(u'that user registers with the email address from the invitation')
def invited_user_registers(context):
    context.active_user = Bot(context.ptmscoutapp)
    context.active_user.email = context.invited_email
    
    context.active_user.register()
    context.active_user.activate()

@when(u'other users search for proteins in my dataset')
def user_search_for_ack1_and_homo_sapiens(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/proteins', status=200)
    form = context.result.form
    
    form.set('acc_search', "ACK1")
    form.set('species',"homo sapiens")
    
    context.result = form.submit()

@when(u'other users lookup proteins that have data in my dataset')
def user_lookup_ack1_homo_sapiens(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/proteins/35546/modifications', status=200)

@when(u'other users access experimental data for proteins that have data in my dataset')
def user_lookup_ack1_homo_sapiens_data(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/proteins/35546/data', status=200)

@when(u'other users attempt to access my experiment directly')
def user_lookup_experiment_26(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/experiments/26')

@when(u'I enter another user email address in "Share dataset"')
@patch('ptmscout.utils.mail.send_automail_message')
def share_experiment_26(context, patch_mail):
    context.owner_user.login()
    context.result = context.ptmscoutapp.get('/account/experiments/26/share', status=200)
    form = context.result.forms[0]
    
    form.set('email', context.active_user.email)
    context.result = form.submit()
    context.result.mustcontain(context.active_user.email)
    
    assert patch_mail.called

@when(u'I press the "publish" button on my experiments page')
def publish_experiment_26(context):
    context.owner_user.login()
    context.result = context.ptmscoutapp.get('/account/experiments', status=200)
    
    forms = context.result.forms__get()
    context.result = context.ptmscoutapp.get(forms['publish26'].action, status=200)
    context.result.mustcontain(strings.publish_experiment_confirm_message)
    
    forms = context.result.forms__get()
    context.result = forms['confirm'].submit()
    context.result.mustcontain(strings.publish_experiment_success_message)
    
    context.owner_user.logout()


@then(u'the unregistered user receives an invitation email')
@patch('ptmscout.utils.mail.send_automail_message')
def check_invite_email_sent(context, patch_mail):
    context.result = context.form.submit()
    context.result.mustcontain(strings.user_invited % context.invited_email)
    
    assertContains(context.invited_email, str(patch_mail.call_args))
    assertContains(context.owner_user.fullname, str(patch_mail.call_args))
    assert patch_mail.called
    

@then(u'I should be able to see the experiment')
def owner_can_see_experiment(context):
    context.owner_user.login()
    
    context.result = context.ptmscoutapp.get('/experiments', status=200)
    context.result.mustcontain("Time-resolved mass spectrometry of tyrosine phosphorylation sites")
    context.result.mustcontain("Effects of HER2 overexpression on cell signaling networks governing proliferation and migration")
    
    context.result = context.ptmscoutapp.get('/experiments/26', status=200)
    context.result.mustcontain("Time-resolved mass spectrometry of tyrosine phosphorylation sites")
    
    context.owner_user.logout()
        
@then(u'other users should not be able to see the experiment')
def other_users_cannot_see_experiment(context):
    context.result = context.ptmscoutapp.get('/experiments', status=200)
    assertDoesNotContain("Time-resolved mass spectrometry of tyrosine phosphorylation sites", context.result.normal_body)
    assertDoesNotContain("Effects of HER2 overexpression on cell signaling networks governing proliferation and migration", context.result.normal_body)
    
    context.result = context.ptmscoutapp.get('/experiments/26', status=200)
    context.result.mustcontain("Forbidden")

@then(u'they should receive a 403 forbidden')
def error_403_forbidden(context):
    context.result.mustcontain("Forbidden")
    
@then(u'my experimental data should not appear in the protein listing')
def ack1_exp_data_not_on_search_page(context):
    p = context.result.pyquery

    table_tag = p('table')
    pytable = helpers.parse_table(p, table_tag)

    assertions.assertEqual(['activated p21cdc42Hs kinase [Homo sapiens].', 'ACK1', 'homo sapiens', '1036', '1', '1', 'Y', 'Phosphotyrosine'], pytable[1])
    assertions.assertEqual(['Activated CDC42 kinase 1; Short=ACK-1; AltName: Full=Tyrosine kinase non-receptor prot', 'TNK2', 'homo sapiens', '1038', '1', '22', 'Y,S,T', 'Phosphoserine, Phosphotyrosine, Phosphothreonine'], pytable[2])
    assertions.assertEqual(len(pytable), 3)


@then(u'my experimental data should not appear in the protein summary')
def ack1_exp_data_not_on_protein_page(context):
    context.result.mustcontain("Y518")
    assertDoesNotContain("Y857", context.result.normal_body)
    assertDoesNotContain("Y858", context.result.normal_body)
    context.result.mustcontain("Quantitative analysis of EGFRvIII cellular signaling networks reveals a combinatorial therapeutic strategy for glioblastoma")
    assertDoesNotContain("Time-resolved mass spectrometry of tyrosine phosphorylation sites", context.result.normal_body)
    assertDoesNotContain("Effects of HER2 overexpression on cell signaling networks governing proliferation and migration", context.result.normal_body)
    
@then(u'that user can now see my specific dataset')
def step_1(context):
    user_lookup_experiment_26(context)
    context.result.mustcontain("Time-resolved mass spectrometry of tyrosine phosphorylation sites")

@then(u'everyone should be able to see my experiment')
def step_2(context):
    context.result = context.ptmscoutapp.get('/experiments/26')
    context.result.mustcontain("Time-resolved mass spectrometry of tyrosine phosphorylation sites")
