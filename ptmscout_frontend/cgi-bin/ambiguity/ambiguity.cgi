#!/usr/local/bin/python
from pathname import *
import cgi
import cgitb
if displayPythonErrors:
	cgitb.enable()
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import sys

sys.path.append(path+'includes')
from template import *
from peptideFunctions import *

from tableFunctions import *
def main():
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    	c=db.cursor()
	form = cgi.FieldStorage()
	default = form.getvalue("default","no")
	if default == "yes": default = True
	else: default = False
	expid=form.getvalue("expid",None)
	printHeader("Report Ambiguity",expid)
	printExpHeader(expid,c)
	peps = getPeps(expid,c)
	ambiguousPeps = getAmbiguousPeps(peps,c,expid)
	msids= list(sets.Set([getMSFromPep(item[0],expid,c) for item in ambiguousPeps]))
	if len(msids)>0:
		print "<a href='%stitle=Ambiguity_Report'>Tutorial for changing protein assignments of ambiguous peptides</a><BR><BR> "%helpPath
		print "<form method='GET' action='newDataFile.txt'>"
		print "<div style='width:100%;text-align:right;display:none;margin-right:2em;' id='export_file'>"
		print "<input type=submit value='Export updated experiment file'>"
		print "</div>"
		makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","ambiguity"],"ambiguity.cgi",msids=alphabetizeMS(msids,c),width='70',icons=False,syns=False,newtab=True,default=default)
		print "<input type=hidden name='expid' value='%s'>"%expid
		print "</form>"
		print "<script>document.getElementById('export_file').style.display='block';</script>"
	else:
		print "<h1>No experiment selected</h1>"
	printFooter()




main()
