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
		Then the user should be sent an email with a link to the dataset which contains:
		  | proteins | peptides | rejected | errors |
		  | 3        | 7        | 1        | 1      |
		And the user should be able to export their dataset with additional annotations:
		  | field                | elements |
		  | nearby_modifications | 17       |
		  | protein_domains      | 37       |
		  | site_regions         | 2        |
		  | site_domains         | 5        |
		  | scansite_kinase      | 15       |
		  | scansite_bind        | 4        |
		  | protein_GO_BP        | 153      |
		  | protein_GO_MF        | 62       |
		  | protein_GO_CC        | 55       |
		And the user should be able to view summary data
		And the user should be able to use the dataset explorer
		And the user should not be able to use other experiment specific tools
		
    Scenario: Export matlab files for MCAM analysis
    	Given a user uploads multiple clusterings for an experiment
    	When the user chooses to export MCAM input files
    	Then the user should be sent an email with a link to download their analysis
    	And the user should be able to download an archive containing the MCAM files
