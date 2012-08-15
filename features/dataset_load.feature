Feature: Dataset Load
	In order to transfer a dataset to PTMScout
	People should be able to load a dataset

	Scenario: Correct Dataset
		Given a correctly formatted dataset of phosphorylation data
		When I press submit dataset
		Then I should be sent an email to a link of my experiment with no errors and contains: 
			17 peptides 
			16 proteins 

	Scenario: More than one pep
		Given a dataset file with more than one peptide column
		When I press submit dataset
		Then I should see an error that I have more than one peptide column and suggests I check for the phrase 'pep' in all columns of my header

	Scenario: No Acc
		Given a dataset with no accession column
		When I press submit dataset
		Then I should see an error that says I have no accession column

	Scenario: Peptide Mismatch	
		Given I submit a dataset that has an incorrect peptide to protein match
		When I press submit dataset 
		Then I should be sent an email to a link of my experiment which contains one error and contains:
			16 peptides
			15 proteins

	Scenario: Bad acc
		Given I submit a dataset that has an accession that looks like a GenPept accession, but is not currently there
		When I submit the dataset
		Then I should be sent an email to a link of my experiment which contains one error and contains:
			17 peptides
			16 proteins
		
	
