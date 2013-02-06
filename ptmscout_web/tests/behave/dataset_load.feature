Feature: Dataset Load
    In order to transfer a dataset accurately to PTMScout
    People should be able to see the errors that occurred during data loading
    
    Scenario: Correct Dataset
          Given a user submits a correctly formatted dataset of phosphorylation data
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 17       | 15       | 0      |
    
    Scenario: Load a correct dataset with degenerate methylation names
          Given a user has loaded a dataset in which modifications have varying specificities of naming
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 7        | 7        | 1      |
    
    Scenario: Load a dataset with bad amino acid assignments
          Given a user has loaded a dataset and the modification type does not match the amino acid for that species
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 7        | 7        | 3      |
          And the experiment browser should contain the correct peptides


    Scenario: Peptide Mismatch    
          Given a user submits a dataset that has an incorrect peptide to protein match
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 16       | 14       | 1      |

    Scenario: Bad acc
          Given a user submits a dataset with an accession that looks like a GenPept accession, but is not a valid accession number
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 17       | 15       | 1      |

    @runme
    Scenario: Multiple modification types in one peptide with assignment
          Given a user submits a dataset and a single peptide has more than one flavor of modification and mod_type assignments match the sequential order of modifcations
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 17       | 15       | 0      |

    Scenario: Multiple modification types in one peptide without assignment
          Given a user submits a dataset in which a single peptide has more than one flavor of modification and the mod_type entry does not match the amino acids
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 16       | 14       | 1      |

    Scenario: Handle a dataset where the same peptide is reported
          Given a user submits a dataset in which the same peptide
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 10       | 7        | 0      |

    Scenario: No peptide passes in dataset
    	Given a user submits a dataset in which no protein accession is found
    	Then the user should be sent an email with a link to the experiment which contains:
    	  | peptides | proteins | errors |
    	  | 0        | 0        | 9      |
    
    Scenario: Handle an isoform
          Given a user submits a dataset in which an isoform specific record is included
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 18       | 16       | 0      |
          
    Scenario: Handle a viral proteins by checking PTM types of host organism
          Given a user submits a dataset in which viral proteins are included
          Then the user should be sent an email with a link to the experiment which contains:
            | peptides | proteins | errors |
            | 6        | 3        | 0      |

