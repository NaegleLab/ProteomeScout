@given(u'I want to create an account and I have a .edu email address')
def step(context):
    assert True

@given(u'that I have registered for an account')
def step(context):
    assert True

@given(u'I have forgotten my password')
def step(context):
    assert True
    
    
@when(u'I enter a my user information')
def step(context):
    assert True

@when(u'I follow the link in my e-mail to activate my account')
def step(context):
    assert True

@when(u'I press "I forgot my password"')
def step(context):
    assert True
    
    
@then(u'I should be sent an e-mail on how to activate my account')
def step(context):
    assert True

@then(u'I should be able to log in to my account')
def step(context):
    assert True

@then(u'I should be sent a temporary token that I can follow to set a new password')
def step(context):
    assert True