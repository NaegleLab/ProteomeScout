# <your-path> location of your clone
# <your-url> your URL base
# REPLACE ALL <xxx> with your specific paths
path = '<your-path>/ptmscout_frontend/cgi-bin/'
graphDir = "<your-url>/cgi-bin/imageScripts/"
urlpath = "<your-url>/cgi-bin/"
entropyPath = "<your-path>/ptmscout_frontend/cgi-bin/tempFiles/entropy/"
tablePath = "<your-path>/ptmscout_frontend/cgi-bin/tempFiles/tables/"
outPath = "<your-path>/ptmscout_frontend/cgi-bin/tempFiles/motifFiles/"
motifPath = "<your-path>/ptmscout_frontend/cgi-bin/subset/motifs/"
clusterPath = "<your-path>/ptmscout_frontend/cgi-bin/tempFiles/clusters/"
urlPath = "<your-url>/cgi-bin/"
summaryPath = "/tmp/ptmscout/"
dataPath = "<your-path>/ptmscout_backend/datasets/"
mpl="<your-path>/ptmscout_frontend/cgi-bin/imageScripts/mpl/"
logPath = "<your-path>/ptmscout_backend/log/"
helpPath = "http://ptmscout.mit.edu/docs/index.php?" 
motifScratch = "/tmp/ptmscout"
imgPath = "<your-url>/images/"
jsPath = "<your-url>"
eggPath = "<your-path>/ptmscout_frontend/cgi-bin/pegg/"
libPath = "<your-path>/ptmscout_backend/ptmscout_source/lib/"
scriptPath = "/data/ptmscout/ptmscout_backend/scripts/dataIO/"
fastaPath = "<your-path>/ptmscout_frontend/cgi-bin/tempFiles/fasta/"
weblogoPath = "<your-path>/ptmscout_frontend/cgi-bin/weblogo/"
dataSummary = "<your-path>/ptmscout_frontend/cgi-bin/tempFiles/datasets/"
scratchPath = "/tmp/ptmscout/"
sysEmail = "your-email@your-institute.edu"

import os
os.environ['PYTHON_EGG_CACHE']=eggPath

import sys
sys.path.append(path+'includes')
from secureDB import *

displayPythonErrors = False

PERL_DB = 'P'; #set this up in ptmscout_backend/ptmscout_source/lib/globalVars.pm
TRACK = 0;
