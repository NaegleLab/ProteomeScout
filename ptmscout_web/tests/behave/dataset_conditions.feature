Feature: Dataset Conditions UI
    Users should be able to describe their datasets using experimental conditions
    Automatic suggestions will come from existing experimental metadata


    Scenario: Describe cell, tissue or condition types
          Given a user is selecting experimental parameters for a dataset
          Then give them options to add cell, tissue, and condition types

    Scenario: Choose condition type
          Given a user is selecting experimental parameters for a dataset
          When they have selected add condition type
          Then give users a selection between drug, stimulation, or environmental conditions

    Scenario: Describe conditions
          Given a user is selecting experimental conditions for a dataset
          Then automatically suggest suitable completions based on existing entries
            | user_input | field_name   | suggestion                     |
            | d          | drug         | dasatinib                      |
            | d          | drug         | doxirubicin                    |
            | H          | cell_type    | HEK 293T                       |
            | H          | cell_type    | HEK 293                        |
            | H          | cell_type    | HELA                           |
            | H          | cell_type    | HEPG2                          |
            | H          | cell_type    | HOP62                          |
            | H          | cell_type    | HOP92                          |
            | H          | cell_type    | HS578T                         |
            | H          | cell_type    | HSG                            |
            | H          | cell_type    | HT1080                         |
            | b          | tissue_type  | BM-CD105+Epithelial cells      |
            | b          | tissue_type  | BM-CD33+Myeloid                |
            | b          | tissue_type  | BM-CD34_                       |
            | b          | tissue_type  | BM-CD70+EarlyErythroid         |
            | b          | tissue_type  | Whole Blood                    |
            | b          | tissue_type  | Whole Brain                    |
            | b          | tissue_type  | Bronchial Epithelial Cells     |
            | b          | tissue_type  | Fetal Brain                    |
            | b          | tissue_type  | Oflactory Bulb                 |
            | b          | tissue_type  | PB-CD19+Bcells                 |
            | b          | tissue_type  | Bone Marrow                    |
            | b          | tissue_type  | Lymphoma Burkitts Daudi        |
            | b          | tissue_type  | Lymphoma Burkitts Raji         |

## BONUS feature -- not high priority ##
#    Scenario: Pre-populate condition fields  
#          Given a user is loading a dataset
#          When they have selected loaded a dataset with data fields that contain condition or cell type data
#          Then pre-populate the add condition with that value 
#            | user_input   | field_name | suggestion |
#            | condition    | EGF        | EGF        |