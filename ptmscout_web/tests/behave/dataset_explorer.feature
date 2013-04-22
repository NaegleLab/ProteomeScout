Feature: Dataset Explorer
    Users should be able to explore uploaded experiments, upload annotations to an experiment, and upload 
    peptides for annotation and exploration

    Scenario: Assign annotations to measured peptides in an experiment
          Given a user uploads annotations to an experiment
          Then the user should be sent an email with a link to the dataset explorer with:
            | annotations | errors   |
            | 201         | 3        |
          And the user should be able to see their annotations in the dataset explorer
          And the user should be able to view the clusterings
          And the user should be able to filter by nominative data
          And the user should be able to filter by numerical data
          And the user should be able to view nominative features in the enrichment data
          
         @runme
	Scenario: Load residues of interest from a dataset and annotate them
		Given a user uploads a file containing non-mass spec experimental data
		Then the user should be sent an email with a link to the experiment which contains:
		  | proteins | peptides | rejected | errors |
		  | 5        | 8        | 1        | 1      |
		And the user should be able to export their dataset with additional annotations:
		  | field                | elements |
		  | modifications        |          |
		  | nearby_modifications |          |
		  | GO_terms             |          |
		  | scansite_kinase      |          |
		  | scansite_bind        |          |
		  | pfam_domains         |          |
		  | pfam_sites           |          |
		And the user should be able to view summary data
		And the user should be able to use the dataset explorer
		And the user should not be able to use other experiment specific tools