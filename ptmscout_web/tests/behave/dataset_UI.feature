Feature: Dataset UI
	In order for users to see and analyze a dataset
	Users should be able to load a dataset 


	Scenario: Load a dataset with uncrecognized headers
		  Given a user is loading a dataset
		  When headers cannot be identified that conform to acc, MOD_TYPE, data, stddev and pep
		  Then show the user their headers and offer drop down menus to assign the necessary headers in the case of data columns the write-in value for type should be propagated to all data columns selected

        Scenario: Loading a dataset with explicit headers
		  Given a user is loading a dataset
		  When data columns exist with declarations acc, MOD_TYPE, pep, data:type:value and matching stddev:type:value
		  Then show the user column assignments and data column assignments with type and value fields and show them an example of the graph that results from that data sequence

	Scenario: Append dataset
		  Given a user is loading a dataset
		  When the dataset exists and new entries are being appended
		  Then pre-populate all fields of data loading 
	   
	   Scenario: Reload a dataset
	   	  Given a user is loading a dataset
		  When the user wants the dataset to replace an alternate dataset		    
		  Then pre-populate all fields of data loading and replace the former experiment 

	Scenario: Load a dataset with bad characters in the peptide column entries
		  Given a user is loading a dataset 
		  When non-alphabetic characters appear in all of the entries
		  Then show the user an error and reset their interface to the initial file loading page 
		  
	Scenario: Load a dataset with modifications described in peptide column
		  Given a user is loading a dataset
		  When all amino acids in lower case do not match any species possibilities for that type of modification then report error and return interface to the intial file loading page

        Scenario: Replicate entries for a modified peptides exist with different data
		  Given a user is loading a dataset with data columns  
		  When identical peptide/protein pairs exist with different data and no explicit run column
		  Show the user that replicate data appears to have been detected and give them an option to select a run column from their header or confirm that replicates can be automatically assigned numerically in order of the appearance

        Scenario: Load a dataset with a pmid 
		  Given a user is loading a datsaet
		  When a pubmed id exists and is given by the user 
		  Automatically fill out the publication information from the pubmed record

	Scenario: Describe experimental conditions
	`	  Given a user is loading a dataset
		  When they are filling out the experimental condition fields
		  Automatically suggest suitable completions based on existing entries where cell lines are constrained by species


	

