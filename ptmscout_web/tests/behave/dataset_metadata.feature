Feature: Dataset metadata
    Users should be able to submit publication information for a dataset 
    
    Scenario: Load a dataset with a pmid
          Given a user is submitting publication information for a new dataset
          When the user submits a pubmed ID to fill out citation information
          Then automatically fill out the publication information from the pubmed record