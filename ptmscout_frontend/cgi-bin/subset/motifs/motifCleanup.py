import time
import os
import sys
import socket
from smtplib import SMTP,SMTPSenderRefused, SMTPRecipientsRefused
from pathname import *
try:
    f = open(motifPath+'numsfile.txt','r')
    lines = f.readlines()
    f.close()
except IOError: lines = []

print "DEBUG: HELLOW WOLRD"

### make sure we have the key
#files = os.listdir('/data/ptmscout/cgi-bin/test/subset/motifs/')
files = os.listdir("%s/subset/motifs/"%path)
if 'motif.lock' not in files:
    key = True
    open('motif.lock','w').write(' ')
else:
    key = False
    print "ERROR, can't open lock file"
    sys.exit(0)
    
####
#w = open('/data/ptmscout/cgi-bin/test/subset/motifs/working.txt','r')
w = open("%s/subset/motifs/numsfile.txt"%path,'r')
working = w.readlines()
workingNums = []
for item in working:
    try:
        workingNums.append(int(item.strip()))
    except ValueError: pass

w.close()
w = open('numsfile.txt','w')

toKeep = []
getRid = []
now = time.time()

##keep if time stamp less than 24 hours prior
for line in lines:
    #print "DEBUG line"
    try:
        try:
            num, then, email = line.split('\t')[0].strip(), line.split('\t')[1].strip(), line.split('\t')[2].strip()
        except IndexError:
            num, then = line.split('\t')[0].strip(), line.split('\t')[1].strip()
            email = ""
        if abs(int(then)-int(now))< 24*60*60:
            toKeep.append((num,then,email))
        else:
            getRid.append(num)
    except ValueError: pass

## for getRids, delete the outHTML, fg, and bg files
for item in getRid:
    #print "Removing item"
    fgf = motifPath+'foreground'+str(item)
    bgf = motifPath+'background'+str(item)
    html = motifPath+'htmls/'+'outHTML'+str(item)
    os.popen("rm -f %s"%fgf)
    os.popen("rm -f %s"%bgf)
    os.popen("rm -f %s"%html)

## for toKeeps, generate new outHTML
## if an email address is given, send the results link to the given address
## if len(toKeep)>2: toKeep = toKeep[0:2]

for item in toKeep:
    #print "DEBUG Keeping item"
    if int(item[0]) not in workingNums:
        if 'outHTML'+str(item[0]) not in os.listdir(motifPath+'htmls/') and len(working)<2:
            w.write(str(item[0])+'\n')
            print "python %ssubset/motifs/makeHTML.py %s"%(path,item[0])
            os.popen("python %ssubset/motifs/makeHTML.py %s"%(path,item[0]))
            if email!= " " and email != "":
                #print "Attempting email"
                x=SMTP("localhost")
                try:
                    x.sendmail(sysEmail,email,"Your PTMScout motif results can be found here: %s."%"%s/subset/motifs/motifFetch.py?num=%s"%(urlpath,str(item[0])))
                except SMTPSenderRefused: pass
                except SMPTRecipientsRefused: pass
            
        
    
####

### rewrite the numfile
f = open(f.name,'w')
for item in toKeep:
    f.write(item[0]+'\t'+item[1]+'\n')
f.close()


#### get rid of all table files
os.system("rm -f %s/subset/motifs/tables/*"%path)

### unlock
os.system('rm -f motif.lock')
w.close()

#### get rid of all table files
os.system("rm -f %s/subset/motifs/tables/*"%path)



