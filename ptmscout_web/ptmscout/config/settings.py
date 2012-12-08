documentationUrl="http://ptmscout.mit.edu/docs/index.php?"
adminEmail = "matt.matlock@gmail.com"
automailerEmail = "automailer@ptmscout.wustl.edu"
pubmedUrl = "www.ncbi.nlm.nih.gov/pubmed/%d"

MINIMUM_PASSWORD_LENGTH = 7


accession_urls = {'swissprot':"http://www.uniprot.org/uniprot/%s",
                  'entrez_protein':"http://www.ncbi.nlm.nih.gov/protein/%s",
                  'gi':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s",
                  'refseq':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s",
                  'GO':"http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s",
                  'genbank':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s"}


ptmscout_scratch_space = "/tmp"
experiment_data_file_path = "data/experiments"
ptmscout_path = "/data/ptmscout/ptmscout_web"

email_regex = "[a-z0-9\.\-\_]+@[a-z0-9\.\-\_]+\.([a-z]+)$"

DISABLE_PFAM = False
DISABLE_SCANSITE = True
DISABLE_QUICKGO = False
DISABLE_PICR = False 
DISABLE_UNIPROT_QUERY = False

valid_domain_suffixes = set(['edu','gov'])






