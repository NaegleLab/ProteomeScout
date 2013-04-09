import os

documentationUrl="http://ptmscout.mit.edu/docs/index.php?"
adminEmail = "ptmscout@seas.wustl.edu"
automailerEmail = "automailer@naegledev.seas.wustl.edu"
pubmedUrl = "www.ncbi.nlm.nih.gov/pubmed/%d"

MINIMUM_PASSWORD_LENGTH = 7


accession_urls = {'swissprot':"http://www.uniprot.org/uniprot/%s",
                  'uniprot':"http://www.uniprot.org/uniprot/%s",
                  'entrez_protein':"http://www.ncbi.nlm.nih.gov/protein/%s",
                  'gi':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s",
                  'refseq':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s",
                  'GO':"http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s",
                  'genbank':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s"}

mod_separator_character = ';'

naegle_lab_url = 'http://naegle.wustl.edu/'
pfam_family_url = 'http://pfam.sanger.ac.uk/family/'
ptmscout_scratch_space = "/tmp"
experiment_data_file_path = "data/experiments"
export_file_path = "data/export"
pfam_map_file_path = "ptmworker/helpers/pfam.map"
ptmscout_path = "/data/ptmscout/ptmscout_web"
motif_script_path = os.path.join('scripts', 'motif')

email_regex = "[a-z0-9\.\-\_]+@[a-z0-9\.\-\_]+\.([a-z]+)$"

DISABLE_PFAM = False
DISABLE_SCANSITE = False
DISABLE_QUICKGO = False
DISABLE_PICR = False
DISABLE_UNIPROT_QUERY = False

JOB_AGE_LIMIT = 7 * 86400
REVIEWER_ACCOUNT_EXPIRATION_DAYS = 365

# rate limits in queries per second
SCANSITE_RATE_LIMIT = 3

valid_domain_suffixes = set(['edu','gov'])


isoform_sequence_diff_pfam_threshold = 50



experiment_files_prefix = "experiment_data"
annotation_files_prefix = "annotation_data"
