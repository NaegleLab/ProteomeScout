Feature: Dataset Conditions UI
    Users should be able to describe their datasets using experimental conditions
    Automatic suggestions will come from existing experimental metadata


    Scenario: Describe cell type
          Given a user is loading a dataset
          When they have selected add cell type 
          Then automatically suggest suitable completions based on existing entries 
            | user_input | field_name | suggestion |
            | H          | cell_type  | HEK 293T, HEK 293, HELA, HEPG2, HOP62, HOP92, HS578T, HSG, HT1080 |
 
    Scenario: Describe tissue type
          Given a user is loading a dataset
          When they have selected add tissue type 
          Then automatically suggest suitable completions based on existing entries 
            | user_input | field_name   | suggestion |
            | b          | tissue_type  | BM-CD105+Epithelial cells, BM-CD33+Myeloid, BM-CD34_, BM-CD70+EarlyErythroid, Whole Blood, Whole Brain, Bronchial Epithelial Cells, Fetal Brain, Oflactory Bulb, PB-CD19+Bcells, Bone Marrow, Lymphoma Burkitts Daudi, Lymphoma Burkitts Raji |

    Scenario: Choose condition type
          Given a user is loading a dataset
          When they have selected add condition type
          Then give users a selection between drug, stimulation, or environmental conditions

    Scenario: Describe condition 
          Given a user is loading a dataset
          When they have selected add condition type equal to drug
          Then automatically suggest suitable completions based on existing entries 
            | user_input | field_name | suggestion |
            | d      | drug       | dasatinib, doxirubicin |

## BONUS feature -- not high priority ##
#    Scenario: Pre-populate condition fields  
#          Given a user is loading a dataset
#          When they have selected loaded a dataset with data fields that contain condition or cell type data
#          Then pre-populate the add condition with that value 
#            | user_input   | field_name | suggestion |
#            | condition    | EGF        | EGF        |