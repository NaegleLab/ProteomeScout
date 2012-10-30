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
          Then show the user that bad peptide strings have been detected
          
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
