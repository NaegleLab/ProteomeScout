Feature: Dataset Conditions UI
    Users should be able to describe their datasets using experimental conditions
    Automatic suggestions will come from existing experimental metadata


    Scenario: Describe cell, tissue or condition types
          Given a user is selecting experimental conditions for a dataset
          Then give the user options to add cell, tissue, drug, stimulation, or environmental conditions

    Scenario: Describe conditions
          Given a user is selecting experimental conditions for a dataset
          Then automatically suggest suitable completions based on existing entries
            | user_input | field_name   | suggestion                     |
            | H          | cell         | HEK 293T                       |
            | H          | cell         | HEK 293                        |
            | H          | cell         | HELA                           |
            | H          | cell         | HEPG2                          |
            | H          | cell         | HOP62                          |
            | H          | cell         | HOP92                          |
            | H          | cell         | HS578T                         |
            | H          | cell         | HSG                            |
            | H          | cell         | HT1080                         |
            | b          | tissue       | BM-CD105+Epithelial cells      |
            | b          | tissue       | BM-CD33+Myeloid                |
            | b          | tissue       | BM-CD34_                       |
            | b          | tissue       | BM-CD70+EarlyErythroid         |
            | b          | tissue       | Whole Blood                    |
            | b          | tissue       | Whole Brain                    |
            | b          | tissue       | Bronchial Epithelial Cells     |
            | b          | tissue       | Fetal Brain                    |
            | b          | tissue       | Oflactory Bulb                 |
            | b          | tissue       | PB-CD19+Bcells                 |
            | b          | tissue       | Bone Marrow                    |
            | b          | tissue       | Lymphoma Burkitts Daudi        |
            | b          | tissue       | Lymphoma Burkitts Raji         |
            | d          | drug         | dasatinib                      |
            | d          | drug         | doxirubicin                    |

## BONUS feature -- not high priority ##
#    Scenario: Pre-populate condition fields  
#          Given a user is loading a dataset
#          When they have selected loaded a dataset with data fields that contain condition or cell type data
#          Then pre-populate the add condition with that value 
#            | user_input   | field_name | suggestion |
#            | condition    | EGF        | EGF        |