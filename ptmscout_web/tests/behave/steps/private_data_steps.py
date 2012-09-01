from behave import *

@given(u'I have logged in with my username and password')
def step(context):
    assert True

@given(u'I want to share my data with another user')
def step(context):
    assert True

@given(u'I want to make my private dataset available to the public')
def step(context):
    assert True



@when(u'I load a dataset and mark it private')
def step(context):
    assert True

@when(u'I enter their username in "Share dataset"')
def step(context):
    assert True

@when(u'I press the "publish" button')
def step(context):
    assert True
    


@then(u'I should be able to see the experiment')
def step(context):
    assert True
    
@then(u'other users should not be able to see the experiment')
def step(context):
    assert True

@then(u'that user can now see my specific dataset')
def step(context):
    assert True

@then(u'everyone should be able to see my experiment')
def step(context):
    assert True