Feature: Error Report
	 In order for a user to easily fix errors that occurred during a dataset load
	 A user should be able given a report and a dataset that only contains the errors
	 
	 Scenario: Public viewing of errors
	 	   Given a user opens an error log
		   When the user is not the owner of the experiment
		   Then the user should see a list of the proteins and peptides that were rejected during data loading

         Scenario: Owner viewing of errors
	 	   Given a user opens an error log
		   When the user is the owner of the experiment
		   Then the user should see a list of the proteins and peptides that were rejected during data loading
		   And the user should be given a link to download the subset of the dataset that did not load with an additional error explanation column
		   And the user should be given a link to append data to the dataset




