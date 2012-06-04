#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import sys

sys.path.append(path+'includes')
sys.path.append(path+'batch')
from template import *
from peptideFunctions import *
from tableFunctions import *

def main():
	
	form = cgi.FieldStorage()
        search = form.getvalue("acc_search","nothing")
	letter=form.getvalue("letter","A")
	stringency = form.getvalue("stringency","Low")
	printHeader("PTMScout, Sample Data Page")
	print "<br><br>"
	exp_id= form.getvalue("exp_id",None)
	db = MySQLdb.connect(user= "mgymrek", passwd="%s"%mysql, db="webdev")
	c=db.cursor()
	
	print '<BR><a href="experiment.cgi?expid=%s">Back to Experiment</a><BR><BR>'%(str(exp_id))

	print """<form name ="search" method="POST" action ="proteinsearch.cgi">
	<input type = "hidden" name = "letter" value = "%s">
	         Search By Protein:&nbsp;&nbsp; <input type="text" name="acc_search"><br><BR>
		 Prediction Stringency: &nbsp;&nbsp; Low <input type="radio" name=stringency value = "Low" checked>&nbsp;&nbsp;
		                                     Medium  <input type="radio" name = stringency value = "Medium">&nbsp;&nbsp;
						     High <input type="radio" name = stringency value = "High">
						     <br><br><input type="submit" value="search">
	</form><BR><BR><BR>""" % (letter)
	printSetOptions(stringency) 
	if search!="nothing":
		makeWholeExpTable(exp_id,c,search,stringency,form,["protein","species","site","pep_aligned"],"sample.cgi",letter=letter,width='75%')
	print '<br><br>'
	db.close()
	print "<BR>"*5
	printFooter()








	
main()

