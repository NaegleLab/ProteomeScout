Feature: Private Data
	Users should be able to login and create accounts
	
	Scenario: Create an account
		Given I want to create an account and I have a .edu email address
		When I enter a my user information 
		Then I should be sent an e-mail on how to activate my account
		
	Scenario: Activate an account
		Given I have registered for an account
		When I follow the link in my e-mail to activate my account
		Then I should be able to log in to my account

	Scenario: Forgot password	
		Given I have an existing account
		When I press "I forgot my password" on the login page and enter my email address
		Then I should be sent a temporary password to login