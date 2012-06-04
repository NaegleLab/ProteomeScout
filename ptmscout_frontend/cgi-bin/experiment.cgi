#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN""http://www.w3.org/TR/html4/loose.dtd">"""
import sys

sys.path.append(path+'includes')
from template import *
def main():
	form = cgi.FieldStorage()
	exp_id = int(form.getvalue('expid',None))
	printHeader("Experiment Home",expid=exp_id)
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
	c=db.cursor()
	query = """select id from experiment where export = 0"""
	c.execute(query)
	x=c.fetchall()
	exps = [str(item[0]) for item in x]
	try:
		printExpHeader(exp_id,c)
		query = """select * from experiment where id = %s"""%exp_id
		c.execute(query)
		x=c.fetchall()
		if len(x)>0:
			print """<fieldset><legend>Tools</legend><dl>"""
			if str(exp_id) in exps:
				print "<b>There is limited functionality for this dataset because export capability has been restricted. This was either a choice of the person who loaded the dataset or it is because this is a data compendium and not an experiment</b>.<BR><BR>"
			query = """select ambiguity from experiment where id = %s"""%exp_id
			c.execute(query)
			x=c.fetchall()
			try: ambiguity = str(x[0][0])
			except IndexError: ambiguity = "0"
			if str(exp_id) not in exps:
				print '<dt><a href="browse/sample.cgi?exp_id=%s">Browse Dataset</a></dt>'%(str(exp_id))
				print "<dd>View sites and proteins in the experiment</dd>"

			print '<dt><a href="summary/summary.cgi?exp_id=%s">Experiment Summary</a></dt>'%(str(exp_id))
			print '<dd>Tables and pie graphs of all annotations for the entire dataset</dd>'
	
			if str(exp_id) not in exps:
				print '<dt><a href="subset/select.cgi?exp_id=%s">Evaluate Subsets</a></dt>'%(str(exp_id))
				print "<dd>Select a foreground from data or metadata, view composition and evaluate enrichment</dd>"
	
				print '<dt><a href="batch/batch.py?exp_id=%s">Compare Datasets</a></dt>'% (str(exp_id))
				print "<dd>View novel sites compared to other datasets</dd>"

			
				if ambiguity == "1":
					print "<dt><a href='ambiguity/ambiguity.cgi?expid=%s'>Report Ambiguity</a></dt>"%(exp_id)
					print """<dd>See all peptide assignments that might be ambiguously assigned.  Explore and change assignments</dd>"""
			
				
			print "</fieldset>"
		else: print "<h3>No experiment selected</h3>"
	except MySQLdb.OperationalError: print "<h3>No experiment selected</h3>"
	printFooter()
	db.close()

main()

