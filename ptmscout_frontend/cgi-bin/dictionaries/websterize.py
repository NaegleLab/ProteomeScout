from pathname import *
import MySQLdb
import pickle
import sys


db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
print "db:%s"%database
c=db.cursor()  
from makeDicts import *
try:
    expid = sys.argv[1]
except IndexError:
    print """Usage: python websterize.py <experiment_id>"""
    sys.exit(0)

print "Making summary dictionary"
d = makeSummaryDict(expid,c)
print "Printing summary dictionary"
f = open('%sdictionaries/savedDicts/dict%s'%(path,expid),'w')
pickle.dump(d,f)
f.close()

print "Making protein dictionary"
d = makeDictFromProteins(getProteins(expid,c),c,stringency="medium")
print "Printing protein dictionary"
f = open('%sdictionaries/savedDicts/prodict%s'%(path,expid),'w')
pickle.dump(d,f)
f.close()

print "Making peptide dictionary"
d = makeDictFromPeps(getPeptides(expid,c),c,stringency="medium")
print "Printing peptide dictionary"
f = open('%sdictionaries/savedDicts/pepdict%s'%(path,expid),'w')
pickle.dump(d,f)
f.close()

d = makeDictFromMSids(getMSids(expid,c),c,stringency="medium")
f = open('%sdictionaries/savedDicts/msdict%s'%(path,expid),'w')
pickle.dump(d,f)
f.close()
