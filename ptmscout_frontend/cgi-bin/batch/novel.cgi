#!/usr/local/bin/python
from pathname import *
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
from batchFcns import *
import sys

sys.path.append(path+'includes')
from template import *
from tableFunctions import *
def main():
    	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    	c=db.cursor()
	form = cgi.FieldStorage()

	
	expid = form.getvalue('expid',None)
	printHeader("Novel Sites", expid)
	published = form.has_key("published")
	if expid == None: print "<h2>No Experiment Selected</h2>"
	else:
		printExpHeader(expid,c)
		
		print """<a name = 'top'></a>"""
		print """<a href="batch.py?exp_id=%s">Back to Comparison Page</a>"""%(expid)
		print """<br><br><A href="%sexperiment.cgi?expid=%s">Back to Experiment Page</a><br><Br>"""%(urlpath,expid)
		if not (published):query = "select id from experiment"
		else: query = "select id from experiment where published=1"
		c.execute(query)
		x=c.fetchall()
		expList = [str(item[0]) for item in x]
		try:
			expList.remove(expid)
		except ValueError: pass

		### get parent, remove
		query = """select experiment_id from experiment where id = %s"""%expid
		c.execute(query)
		x=c.fetchall()
		if len(x)>0:
			parent = str(x[0][0])
			if parent in expList:
				expList.remove(parent)
		else: parent = None
		### get others that share that parent, remove
		if parent != None and int(parent)!=0:
			query = """select id from experiment where experiment_id = %s"""% parent
			c.execute(query)
			x=c.fetchall()
			if len(x)>0:
				x=[str(item[0]) for item in x]
				for item in x:
					if item in expList: expList.remove(item)

		### get any children, remove
		query = """select id from experiment where experiment_id = %s"""%expid
		c.execute(query)
		x=c.fetchall()
		if len(x)>0:
			x=[str(item[0]) for item in x]
			for item in x:
				if item in expList:
					expList.remove(item)
		exp1Peps = Set(getPeps(expid,c))
		ambiguousPeps=Set([item[0] for item in getAmbiguousPeps(exp1Peps,c,expid)])
		for ex in expList:
			
			exp2Peps = Set(getPeps(ex,c))
			ambiguousPeps=ambiguousPeps.difference(exp2Peps)
			exp1Peps = exp1Peps.difference(exp2Peps)
			amb = getAmbiguousPeps(exp1Peps,c,expid)
        		for item in amb:
            			
            			possible = getPossiblePeps(item[1],c)
				if Set(possible).intersection(Set(exp2Peps))!=None:
            		
					
					if item[0] in exp1Peps:
						## if ANY possible ambiguity resolution options in other dataset, remove it from novel sites
						exp1Peps.remove(item[0])
						ambiguousPeps.add(item[0])
		
		
        	print "<h3>Novel Sites (%s)</h3>"%len(exp1Peps)
		msids = getMSids_batch(exp1Peps,c,expid)
		newms = []
		for item in msids:
			acc = getTableValues(item, c)[1]
			newms.append(item)
		print '<table width="60%%"><tr><td>'
		
		print '</center><BR>'
		if len(newms)>0:
			makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","pep_aligned"],"",msids=alphabetizeMS(newms,c),newtab=True)
			print '<BR><BR><BR>'
		
		print '</td></tr></table>'

		if len(ambiguousPeps)>0:
			print "<h3>Ambiguous Sites (%s)</h3>"%len(ambiguousPeps)
			print "The following sites have ambiguous protein assignments, some of which may indicate a novel site. <a href='%sambiguity/ambiguity.cgi?expid=%s'><br>Change assignments of ambiguous peptides.</a>"%(urlpath,expid)
			msids = getMSids_batch(ambiguousPeps,c,expid)
			newms = []
			for item in msids:
				acc = getTableValues(item, c)[1]
				newms.append(item)
		
			print '<table width="60%%"><tr><td>'
			print '</center><BR>'
			if len(newms)>0:
				makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","pep_aligned"],"",msids=alphabetizeMS(newms,c),icons=False,syns=False,newtab=True)
				
			print '<BR><BR><BR>'
		
			print '</td></tr></table>'
    	printFooter()

	

main()
