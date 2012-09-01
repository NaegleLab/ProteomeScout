Feature: Private Data
	In order for users to analyze data prior to publication 
	Users should be able to load a dataset that only specified credentials can see

	Scenario: Load private data
		Given I have logged in with my username and password
		When I load a dataset and mark it private
		Then I should be able to see the experiment
		And other users should not be able to see the experiment

	Scenario: Share data
		Given I want to share my data with another user
		When I enter their username in "Share dataset"
		Then that user can now see my specific dataset

	Scenario: Publish data
		Given I want to make my private dataset available to the public
		When I press the "publish" button 
		Then everyone should be able to see my experiment