homeUrl="../cgi-bin"
documentationUrl="http://ptmscout.mit.edu/docs/index.php?"
adminEmail = "matt.matlock@gmail.com"
automailerEmail = "automailer@ptmscout.wustl.edu"
pubmedUrl = "www.ncbi.nlm.nih.gov/pubmed/%d"

MINIMUM_PASSWORD_LENGTH = 7

experiment_data_file_path = "data/experiments"

accession_urls = {'swissprot':"http://ca.expasy.org/uniprot/%s",
                  'entrez_protein':"http://www.ncbi.nlm.nih.gov/protein/%s",
                  'gi':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s",
                  'refseq':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s",
                  'GO':"http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s",
                  'genbank':"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=%s"}


ptmscout_path = "/data/ptmscout/ptmscout_web"

email_regex = "[a-z0-9\.\-\_]+@[a-z0-9\.\-\_]+\.([a-z]+)$"

DISABLE_PFAM = False
DISABLE_SCANSITE = False
DISABLE_QUICKGO = False
DISABLE_PICR = False 

valid_domain_suffixes = set(['edu','gov'])


