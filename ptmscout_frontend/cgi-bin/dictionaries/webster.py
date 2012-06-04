from pathname import *
import MySQLdb
import pickle
import sys

db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
c=db.cursor()  
from makeDicts import *
try:
    expid = sys.argv[1]
except IndexError:
    print """Usage: python webster.py <experiment_id>"""
    sys.exit(0)
d = makeSummaryDict(expid,c)
f = open('%sdictionaries/savedDicts/dict%s'%(path,expid),'w')
pickle.dump(d,f)
f.close()
