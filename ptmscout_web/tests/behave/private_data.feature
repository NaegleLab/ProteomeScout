Feature: Private Data
	In order for users to analyze data prior to publication 
	People should be able to login and load a dataset that only their credentials can see

	Scenario: Create an account
		Given I want to create an account and I have a .edu email address
		When I enter a user name, e-mail and password 
		Then I should be given an account

	Scenario: Forgot password	
		Given I have forgotten my password
		When I press "I forgot my password"
		Then I should be sent a temporary token that I can follow to set a new password 

	Scenario: Load private data
		Given I have logged in with my username and password
		When I load a dataset and mark it private	
		Then only I should be able to see the experiment and my experiment should not be listed in the evidence table for data in my experiment

	Scenario: Share data
		Given that I want to share my data with a group member or a reviewer
		When I enter that users name in "Share dataset"
		Then that user can now see my specific dataset

	Scenario: Publish data
		Given that I want to make my private dataset available to the public
		When I press the "publish" button 
		Then everyone should be able to see my experiment and it will show up in the evidences table
