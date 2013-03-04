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
          And the user should be able to filter by numerical data
          And the user should be able to filter by nominative data