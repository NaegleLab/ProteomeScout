#!/usr/local/bin/python
from pathname import *
import cgi
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import MySQLdb
mysql="melissa_mysql"
db=MySQLdb.connect(user="mgymrek",passwd="%s"%mysql,db="webdev")
c=db.cursor()
def getRealParent(expid,c):
	query = """select experiment_id from experiment where id = %s"""%expid
	c.execute(query)
	try:
		x=c.fetchall()[0][0]
		if int(x)==0: return str(expid)
		else: return str(x)
	except IndexError: return str(expid)

import os

def returnAccType(accession):
	f=os.popen("""perl -I"/programs/lib/knaegle" entrez.pl '%s'"""%accession)
	return f.read()
print 'here'
print returnAccType("gi|1234134134")
