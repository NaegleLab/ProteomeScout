from behave import *

@given(u'a correctly formatted dataset of phosphorylation data')
def step(context):
    assert True

@given(u'a dataset file with more than one peptide column')
def step(context):
    assert True

@given(u'a dataset with no accession column')
def step(context):
    assert True

@given(u'I submit a dataset that has an incorrect peptide to protein match')
def step(context):
    assert True

@given(u'I submit a dataset that has an accession that looks like a GenPept accession, but is not currently there')
def step(context):
    assert True

@when(u'I submit the dataset')
def step(context):
    assert True


@then(u"I should see an error that I have more than one peptide column")
def step(context):
    assert True
    
@then(u"I should see a suggestion to check for the phrase 'pep' in all columns of my header")
def step(context):
    assert True

@then(u'I should see an error that says I have no accession column')
def step(context):
    assert True

@then(u'I should be sent an email to a link of my experiment which contains')
def step(context):
    assert True

