#!/usr/local/bin/python
from pathname import *
import cgi
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import cStringIO
import MySQLdb
import os
import sys

sys.path.append(path+'includes')
sys.path.append(path)
import cgitb
if displayPythonErrors:
	cgitb.enable()
from template import printHeader, printFooter
import time
import copy
from sets import Set
from batchFcns import *
def main():
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    	c=db.cursor()
	form = cgi.FieldStorage()
	exp1 = form.getvalue("exp1",None)
	exp2 = form.getvalue("exp2",None)
	expString = form.getvalue('expString',None)
	common = Set()
	if exp1 !=None and exp2 != None:
		common = Set(getPeps(exp1,c)).intersection(Set(getPeps(exp2,c)))
	elif expString != None:
		exps = [int(item) for item in expString.split(':')]
		for item in exps:
			if len(common)==0: common = Set(getPeps(item,c))
			else:
				common = common.intersection(Set(getPeps(item,c)))
		
	print "<html><head><title>PTMScout - Site Comparison</title></head><body>"
	print "<h2>Site Comparison</h2>"
	### print out the two experiments being compared
	if exp1 != None and exp2 != None:
		query = """select name from experiment where id = %s or id = %s"""%(exp1, exp2)
		c.execute(query)
		x=c.fetchall()
		names = [item[0] for item in x]
		for item in names:
			print "Experiment: " +item +"<br>"
		print "<br>"
		print "<BR>Sites common to both experiments:<br><BR>"
	### print out p-sites in common:
	for pep in common:
		query = """select site_type, site_pos ,acc_gene from phosphopep join protein on phosphopep.protein_id=protein.id where phosphopep.id = %s"""%pep
		c.execute(query)
		x=c.fetchall()
		print x[0][2]+' - '+str(x[0][0])+str(x[0][1])+"<br>" 
	print "</body></html>"
	db.close()

main()
