Feature: Dataset Peptide Ambiguity UI
    Users should be able to reassign their peptides using protein accessions from Uniprot

    Scenario: Assign most populated protein accessions to peptides from ambiguous dataset
          Given a user chooses to assign default accessions for an experiment
          When the user submits the extension experiment for loading
          Then the user should be sent an email with a link to the extension experiment which contains:
            | peptides | proteins | rejected | errors |
            | 61       | 45       | 7        | 28     |

