#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""

import sys
sys.path.append(path+'includes')
from template import *
from proteinFcns import *
def main():
	print "<html>"

	form = cgi.FieldStorage()
	db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
	c=db.cursor()
	groupings = getExpGroupings(c)
	pid = form.getvalue("pid",None)
	printHeader("Experimental Measurements", protein_id=pid)
	if pid == None:
		print "No protein selected"
		printFooter()
		sys.exit(0)
	prodict = makeProteinInfoDict(pid,c)
	makeProteinDetails(prodict)
	makeGeneOntologies(prodict,c)
	### get experiments we want, the intersection of the two
	query = """select name,experiment.id from experiment join MS on MS.experiment_id=experiment.id where export = 1 and protein_id=%s group by experiment.id"""%pid
	c.execute(query)
	x=c.fetchall()
	exps = [(item[1],item[0]) for item in x]
#	exps = [item for item in exps if item[0] in groupings.keys()]
	###
	if len(exps)==0: print "No data measurements exist for this protein"
	for item in exps:
		query="""select * from data join MS on data.MS_id=MS.id where experiment_id=%s and protein_id=%s"""%(str(item[0]),str(pid))
		c.execute(query)
		x=c.fetchall()
		if len(x)>0:
			print "<br><BR><BR>"
			print "<a href='%sexperiment.cgi?expid=%s'>%s</a>"%(urlpath,item[0],item[1])
			print "<BR>"
			print "<img alt='graph' height='600' src='%smakeGraph.py?exp_id=%s&pid=%s'>"%(graphDir,item[0],pid)
			print "<BR>"

	print "</html>"

	printFooter()




main()
