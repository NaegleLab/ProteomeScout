Feature: Private Data
	In order for users to analyze data prior to publication 
	Users should be able to load a dataset that only specified credentials can see

	Scenario: Load private data
		Given I have logged in with my username and password
		When I load a dataset and mark it private
		Then I should be able to see the experiment
		And other users should not be able to see the experiment
	
	Scenario: Private experiments should not be directly accessible
		Given I have loaded a dataset and marked it private
		When other users attempt to access my experiment directly
		Then they should receive a 403 forbidden
		
	Scenario: Private data does not appear in protein search
		Given I have loaded a dataset and marked it private
		When other users search for proteins in my dataset
		Then my experimental data should not appear in the protein listing
		
	Scenario: Private data does not appear on protein summary pages
		Given I have loaded a dataset and marked it private
		When other users lookup proteins that have data in my dataset
		Then my experimental data should not appear in the protein summary

	Scenario: Share data
		Given I have loaded a dataset and marked it private
		When I enter another user email address in "Share dataset"
		Then that user can now see my specific dataset

	Scenario: Publish data
		Given I have loaded a dataset and marked it private
		When I press the "publish" button on my experiments page
		Then everyone should be able to see my experiment