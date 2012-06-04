import os
from pathname import *

os.system("rm -f %s/subset/motifs/foreground*"%path)
os.system("rm -f %s/subset/motifs/background*"%path)
os.system("rm -f %s/subset/motifs/log.foreground*"%path)
os.system("rm -f %s/subset/motifs/weboutput*"%path)
os.system("rm -f %s/subset/motifs/working.txt"%path)
os.system("touch %s/subset/motifs/working.txt"%path)

os.system("rm -f %s/tempFiles/clusters/*"%path)
os.system("rm -f %s/tempFiles/entropy/*"%path)
os.system("rm -f %s/tempFiles/fasta/*"%path)
os.system("rm -f %s/tempFiles/motifFiles/*"%path)
os.system("rm -f %s/tempFiles/tables/*"%path)

os.system("rm -f %s/*"%motifScratch)
