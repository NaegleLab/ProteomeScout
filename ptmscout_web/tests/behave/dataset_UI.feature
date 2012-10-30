Feature: Dataset UI
    In order for users to see and analyze a dataset
    Users should be able to load a dataset 


    Scenario: Loading a dataset with unrecognized headers
          Given a user is loading a dataset
          When headers cannot be identified that conform to acc, MOD_TYPE, data, stddev and pep
          Then show the user their headers 
#          And offer drop down menus to assign the necessary headers

    Scenario: Loading a dataset with explicit headers
          Given a user is loading a dataset
          When data columns exist with declarations acc, MOD_TYPE, pep, data:type:value and matching stddev:type:value
          Then show the user column assignments and data column assignments with type and value fields
#          And show the user an example of the graph from their data assignments

    Scenario: Append to a dataset
          Given a user is loading a dataset
          When the dataset exists and new entries are being appended
          Then pre-populate all fields of data loading with assignments from original dataset
       
    Scenario: Reload a dataset
          Given a user is loading a dataset
          When the user wants the dataset to replace another dataset
          Then pre-populate all fields of data loading with assignments from original dataset
          
    Scenario: Extend a dataset
          Given a user is loading a dataset
          When the user wants the dataset to extend another dataset
          Then pre-populate all fields of data loading with assignments from original dataset

    Scenario: Load a dataset with bad characters in the peptide column entries
          Given a user is loading a dataset
          When non-alphabetic characters appear in some of the entries
          Then report the errors and request user confirmation to continue or return to data assignment
          
    Scenario: Load a dataset with modifications described in the peptide column
          Given a user is loading a dataset
          When some amino acids in lower case do not match any species possibilities for that type of modification
          Then show the user that incorrect modifications types have been detected
#          And show the user an option to continue or return to data assignment          

    Scenario: Replicate entries for modified peptides exist with different data
          Given a user is loading a dataset  
          When identical peptide/protein pairs exist with different data and no explicit run column
          Then show the user that replicate data appears to have been detected
#          And show the user an option to select a run column from their header or confirm that replicates can be automatically assigned numerically in order of the appearance

    Scenario: Load a dataset with a pmid
          Given a user is loading a dataset
          When the user submits a pubmed ID to fill out citation information
          Then automatically fill out the publication information from the pubmed record

## The below scenarios relate to the experimental metadata.  Automatic suggestions will come from existing experimental metadata and so far I am adding autocomplete suggestions based on gene expression tables which we can use to pre-populate a metadata table. 
#    Scenario: Describe cell type
#          Given a user is loading a dataset
#          When they have selected add cell type 
#          Then automatically suggest suitable completions based on existing entries 
#	   | user_input | field_name | suggestion |
#	   | H		| cell_type  | HEK 293T, HEK 293, HELA, HEPG2, HOP62, HOP92, HS578T, HSG, HT1080 |
 
#    Scenario: Describe tissue type
#          Given a user is loading a dataset
#          When they have selected add tissue type 
#          Then automatically suggest suitable completions based on existing entries 
#	   | user_input | field_name | suggestion |
# 	   | b		| tissue_type  | BM-CD105+Epithelial cells, BM-CD33+Myeloid, BM-CD34_, BM-CD70+EarlyErythroid, Whole Blood, Whole Brain, Bronchial Epithelial Cells, Fetal Brain, Oflactory Bulb, PB-CD19+Bcells, Bone Marrow, Lymphoma Burkitts Daudi, Lymphoma Burkitts Raji |

#    Scenario: Choose condition type
#          Given a user is loading a dataset
#          When they have selected add condition type
#          Then give users a selection between drug, stimulation, or environmental conditions

#    Scenario: Describe condition 
#          Given a user is loading a dataset
#          When they have selected add condition type equal to drug
#          Then automatically suggest suitable completions based on existing entries 
#	   | user_input | field_name | suggestion |
#          | d		| drug	     | dasatinib, doxirubicin |

## BONUS feature -- not high priority ##
#    Scenario: Pre-populate condition fields  
#          Given a user is loading a dataset
#          When they have selected loaded a dataset with data fields that contain condition or cell type data
#          Then pre-popolate the add condition with that value 
#	   | user_input | field_name | suggestion |
#          | condition		| EGF	     | EGF |