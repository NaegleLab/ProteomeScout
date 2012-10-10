Feature: Dataset Load
	In order to transfer a dataset to PTMScout
	People should be able to load a dataset

	Scenario: Correct Dataset
		Given a user submits a correctly formatted dataset of phosphorylation data
        Then the user should be sent an email with a link to the experiment which contains:
			| peptides | proteins | errors |
			| 17       | 16       | 0      |

	Scenario: More than one pep
		Given a dataset file with more than one peptide column
		When I submit the dataset
		Then the user should see an error that says "data set has more than one peptide column"
		And the user should see the text "check for the phrase 'pep' in all columns of the experiment header"

	Scenario: No Acc
		Given a user submits a dataset with no accession column
		Then the user should see an error that says "data set has no accession column"

	Scenario: Peptide Mismatch	
		Given a user submits a dataset that has an incorrect peptide to protein match
		Then the user should be sent an email with a link to the experiment which contains:
			| peptides | proteins | errors |
			| 17       | 16       | 1      |

	Scenario: Bad acc
		Given a user submits a dataset that has an accession that looks like a GenPept accession, but is not currently there
		Then the user should be sent an email with a link to the experiment which contains:
			| peptides | proteins | errors |
			| 17       | 16       | 1      |
		
	
