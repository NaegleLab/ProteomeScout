from pathname import *
import MySQLdb
import pickle
import sys


db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
print "db:%s"%database
c=db.cursor()
cutoff=500
from makeDicts import *
query = """select experiment.id, count(*) as NUM from experiment join MS on MS.experiment_id=experiment.id where export=1 group by MS.experiment_id having NUM>=%s"""%cutoff
c.execute(query)
x=c.fetchall()
exps = [item[0] for item in x]
for exp in exps:
    print "Experiment %s"%exp
    os.system("python websterize.py %s"%exp)
#d = makeSummaryDict(expid,c)
#f = open('%sdictionaries/savedDicts/dict%s'%(path,expid),'w')
#pickle.dump(d,f)
#f.close()

#d = makeDictFromProteins(getProteins(expid,c),c,stringency="medium")
#f = open('%sdictionaries/savedDicts/prodict%s'%(path,expid),'w')
#pickle.dump(d,f)
#f.close()

#d = makeDictFromPeps(getPeptides(expid,c),c,stringency="medium")
#f = open('%sdictionaries/savedDicts/pepdict%s'%(path,expid),'w')
#pickle.dump(d,f)
#f.close()

#d = makeDictFromMSids(getMSids(expid,c),c,stringency="medium")
#f = open('%sdictionaries/savedDicts/msdict%s'%(path,expid),'w')
#pickle.dump(d,f)
#f.close()
