Feature: Error Report
	 In order for a user to easily fix errors that occurred during a dataset load
	 A user should be able given a report and a dataset that only contains the errors
	 
	 Scenario: Public viewing of errors
	 	   Given a user opens an error log
		   When the user is not the owner of the experiment
		   Then the public user should see a list of the proteins and peptides that were rejected during data loading, if any

         Scenario: Owner viewing of errors
	 	   Given a user opens an error log
		   When the user is the owner of the experiment
		   Then the user should see the public list and be given a link to download the subset of the dataset that did not load, with helpful hints added for each line, and be given a link to immediately append data to that dataset




