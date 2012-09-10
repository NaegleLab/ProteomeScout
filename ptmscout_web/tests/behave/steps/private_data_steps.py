from behave import *
from bot import Bot
from assertions import assertDoesNotContain, assertContains

@given(u'I have logged in with my username and password')
def login_the_active_user(context):
    context.active_user.login()

@given(u'I have loaded a dataset and marked it private')
def create_owning_user_with_datasets(context):
    owner_user = Bot(context.ptmscoutapp)
    owner_user.register()
    owner_user.activate()
    exp_ids=[26,28]
    owner_user.acquire_experiments(exp_ids)
    
    context.owner_user = owner_user




@when(u'I load a dataset and mark it private')
def load_private_dataset(context):
    assert True
    
@when(u'other users search for proteins in my dataset')
def user_search_for_ack1_and_homo_sapiens(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/protein/search', {'protein':"ACK1", 'stringency':1, 'species':"homo sapiens"}, status=200)

@when(u'other users lookup proteins that have data in my dataset')
def user_lookup_ack1_homo_sapiens(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/protein/35546?species=homo%20sapiens', status=200)

@when(u'other users attempt to access my experiment directly')
def user_lookup_experiment_26(context):
    context.active_user.login()
    context.result = context.ptmscoutapp.get('/experiments/26')

@when(u'I enter another user email address in "Share dataset"')
def share_experiment_26(context):
    context.owner_user.login()
    context.result = context.ptmscoutapp.get('/account/experiments/26/share', status=200)
    form = context.form
    
    form.set('email', context.active_user.email)
    context.result = form.submit()
    context.result.mustcontain(context.active_user.email)

@when(u'I press the "publish" button on my experiments page')
def publish_experiment_26(context):
    context.owner_user.login()
    context.result = context.ptmscoutapp.get('/account/experiments/26', status=200)
    forms = context.result.forms__get()
    
    context.result = forms['publish26'].submit()
    context.result.form.submit()
    context.result.mustcontain("You have successfully published this experiment.")




@then(u'I should be able to see the experiment')
def owner_can_see_experiment(context):
    assert True
    
@then(u'other users should not be able to see the experiment')
def other_users_cannot_see_experiment(context):
    assert True

@then(u'they should receive a 403 forbidden')
def error_403_forbidden(context):
    context.result.mustcontain("Forbidden")
    
@then(u'my experimental data should not appear in the protein listing')
def ack1_exp_data_not_on_search_page(context):
    assertContains("Y284", context.result.normal_body)
    assertContains("Y518", context.result.normal_body)
    assertDoesNotContain("Y857", context.result.normal_body)
    assertDoesNotContain("Y858", context.result.normal_body)
    
@then(u'my experimental data should not appear in the protein summary')
def ack1_exp_data_not_on_protein_page(context):
    assertContains("Y284", context.result.normal_body)
    assertContains("Y518", context.result.normal_body)
    assertDoesNotContain("Y857", context.result.normal_body)
    assertDoesNotContain("Y858", context.result.normal_body)
    assertContains("Quantitative analysis of EGFRvIII cellular signaling networks reveals a combinatorial therapeutic strategy for glioblastoma", context.result.normal_body)
    assertDoesNotContain("Time-resolved mass spectrometry of tyrosine phosphorylation sites", context.result.normal_body)
    assertDoesNotContain("Effects of HER2 overexpression on cell signaling networks governing proliferation and migration", context.result.normal_body)
    
@then(u'that user can now see my specific dataset')
def step_1(context):
    user_lookup_experiment_26(context)

@then(u'everyone should be able to see my experiment')
def step_2(context):
    user_lookup_experiment_26(context)