email_regex = "[a-z0-9\.\-\_]+@[a-z0-9\.\-\_]+\.([a-z]+)$"

accession_type_strings = {'ipi':"International Protein Index", 
                     'gene_synonym':"Gene Synonym",
                     'refseq': "RefSeq",
                     'swissprot': "SwissProt",
                     'genbank': "GenBank",
                     'gi': "Gene Index"}

success_header = "Success"

change_password_page_title = "Change Password"
change_password_success_message = "Password successfully changed."

about_page_title = "About PTMScout"

account_activation_page_title = "Account Activation"
account_activation_success_header = "Account Activation Succeeded"
account_activation_failed_header = "Account Activation Failure"
account_activation_success_message = "Your account is now active. Please <a href=\"%s\">login</a>"
account_activation_failed_message = "The specified account is not valid, please try <a href=\"%s\">registering</a>"

account_management_page_title = "Account Management"

experiments_page_title = "Home - Experiments"
experiment_page_title = "Experiment: %s"
experiment_summary_page_title = "Experiment Summary: %s"
experiment_browse_page_title = "Browsing Experiment: %s"
experiment_prediction_page_title = "Experiment Predictions: %s"
experiment_pfam_page_title = "Experiment Protein Families: %s"

experiment_upload_finished_subject = "PTMScout experiment upload completed" 
experiment_upload_finished_message = \
"""PTMScout has finished processing the upload of your experiment: '%s'

Upload Results:

Peptides: %d
Proteins: %d
Errors: %d

You may view the error log for this upload at %s.

Thanks for using PTMScout,
-The PTMScout Team"""
experiment_upload_configure_page_title = "Configure Experiment Data"
experiment_upload_configure_message = "Please verify and assign designations for the column types present in your dataset. Some column types may have been automatically inferred from your headers."

experiment_upload_confirm_page_title = "Confirm Submission"
experiment_upload_confirm_message = "Experiment headers have parsed successfully. You may review warnings below. Are you sure you wish to proceed with the upload?"

experiment_upload_started_page_title = "Upload Started"
experiment_upload_started_message = \
"""Experiment upload process successfully started. An e-mail will be sent to you when the experiment upload is complete. Additionally, you may check the status of the upload by visiting <a href=\"%s\">this page</a>"""


forgotten_password_page_title = "Forgotten Password Retrieval"
forgotten_password_success_header = "Password Reset Success"
forgotten_password_success_message = "Your username and a temporary password have been sent to your e-mail address"
forgotten_password_email_subject = "PTMScout password reset"
forgotten_password_email_message = \
"""%s,

Your password in PTMScout has been reset, your new login credentials are:
Username: %s
Password: %s

Please visit <a href="%s">PTMScout</a> to login.
After logging in, your can change your password <a href="%s">here</a>.

-PTMScout Administrator
"""

login_page_title = "Login"
login_page_success_header = "Login Successful"
login_page_success_message = "You have successfully logged in."

logout_page_title = "Logout"
logout_page_header = "Logout Successful"
logout_page_message = "You have successfully logged out."

my_experiments_page_title = "My Experiments"
publish_experiment_page_title = "Publish Experiment"
publish_experiment_confirm_message = "Are you sure you want to publish this data?"
publish_experiment_success_message = "You have successfully published this experiment."
publish_experiment_already_message = "Experiment has already been published."

privatize_experiment_page_title = "Set Experiment Private"
privatize_experiment_confirm_message = "Are you sure you want to make this data private?"
privatize_experiment_success_message = "You have successfully made this experiment private."
privatize_experiment_already_message = "Experiment has already been made private."

protein_data_page_title = "Experiment Measurements"
protein_ontology_page_title = "Gene Ontologies"
protein_expression_page_title = "Protein Expression"
protein_modification_sites_page_title = "Protein Modification Sites"
protein_search_page_title = "Protein Search"

share_experiment_page_title = "Share Experiment"

upload_page_title = "Upload"
upload_page_header = "Data Upload"

user_invite_page_title = "Invite User"
user_invite_confirm = "User %s is not a registered user, are you sure you wish to invite this user?"
user_invite_email_required = "Email address is required"
user_invited = "An invitation to view your dataset has been sent to %s."
user_invite_email_subject = "PTMScout user %s has invited you to share a dataset"
user_invite_email_message = """
%s,

User %s has invited you to view their dataset '%s', available through PTMScout.

Please <a href=\"%s\">visit</a> to access and view this data.

Thanks,
-The PTMScout Team
"""


user_registration_page_title = "User Registration"
user_registration_success_header = "Registration Successful"
user_registration_success_message = "A confirmation e-mail has been sent to the specified e-mail address. Please check your e-mail to complete your registration."

user_registration_email_subject = "PTMScout Account Activiation Details"
user_registration_email_message = """
%s,

Thank you for choosing PTMScout for your research.

You can activate your new account by visiting <a href=\"%s/activate_account?username=%s&token=%s\">this link</a>.

Thanks,
-The PTMScout Team
"""

view_page_title = "PTMScout Terms of Use"

failure_reason_experiment_file_not_enough_columns = "Not enough columns detected in data file. Verify that the file is TSV format and try again"
failure_reason_experiment_header_no_peptide_column = "Data file did not contain a peptide column"
failure_reason_experiment_header_multiple_peptide_column = "Data file contained multiple peptide columns. Check for the phrase 'pep' in all columns of the experiment header."
failure_reason_experiment_header_no_acc_column = "Data file did not contain an accession column"
failure_reason_experiment_header_multiple_acc_columns = "Data file contained multiple accession columns"

failure_reason_terms_of_use_not_accepted = "You must agree to the terms of use at the bottom of this form before submitting a dataset"
failure_reason_field_value_not_valid = "Field '%s' has invalid value"
failure_reason_field_must_be_numeric = "Field '%s' must be numeric"
failure_reason_required_fields_cannot_be_empty = "Required form field '%s' cannot be empty"
 
failure_reason_form_fields_cannot_be_empty = "Form fields cannot be empty"
failure_reason_new_passwords_not_matching = "Password confirmation did not match"
failure_reason_incorrect_password = "Supplied password was incorrect"

failure_reason_inactive_account = "Account has not been activated"
failure_reason_incorrect_credentials = "Credentials incorrect"
failure_reason_username_inuse = "Username is already in use"
failure_reason_email_not_valid = "Email address is invalid"
failure_reason_email_not_academic = "Email address must belong to .edu domain"
failure_reason_password_too_short = "Password must be at least %d characters in length"
failure_reason_email_address_not_on_record = "E-mail address does not match any user record"


error_resource_not_found_page_title = "Resource Error"
error_protein_not_found_message = "No protein resource exists with the specified ID"

prediction_type_map = {'scansite': "Scansite",
                       'scansite_bind': "Scansite Bind",
                       'scansite_kinase': "Scansite Kinase"}



error_resource_not_ready_page_title = "Resource not ready"
error_resource_not_ready_message = "The resource you are trying to access is currently in processing. If you are the creator of this resource, please visit <a href=\"%s\">this page</a> to check on its status."



 


