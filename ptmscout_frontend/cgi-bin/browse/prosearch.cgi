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
import cgitb
if displayPythonErrors:
	cgitb.enable()
def main():
	
	form = cgi.FieldStorage()
        search = form.getvalue("acc_search","nothing")
	letter=form.getvalue("letter","A")
	stringency = form.getvalue("stringency","Low")
	printHeader("Protein Search")
	print "<br><br>"
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
	c=db.cursor()
	

	print """<form name ="search" method="POST" action ="prosearch.cgi">
	<fieldset><legend>Search</legend>
	<input type = "hidden" name = "letter" value = "%s">
	<p>
	<label for="acc_search">Protein:</label>""" %(letter)
	textboxForSearch("acc_search", form);
	print """
        </p>
	<p>
	<label for="stringency">Prediction Stringency:</label> <br />"""
	
	radioForSearch("stringency", "Low", form)
	print "Low <br />"

	radioForSearch("stringency", "Medium", form)
	print "Medium <br />"

	radioForSearch("stringency", "High", form)
	print "High <br />"
	print """
	<p>
	<label for="species">Species:</label> <br />"""
	selectBoxForSearch("species",["all"]+getAllSpecies(c),form)

	print """
        </p>
	<input type="submit" value="Search">
	</fieldset>
	</form>"""
	printSetOptions(stringency)
	if search!="nothing" or form.getvalue("pname",None) !=None:
		makeWholeExpTable(None,c,search,stringency,form,["MSid","sequence","protein","species","site","pep_aligned","predictions"],"sample.cgi",letter=letter,width='75%',speciesSearch = form.getvalue("species","all"),exactProName=form.getvalue("pname",None))
	print '<br><br>'
	db.close()
	print "<BR>"*5
	printFooter()








	
main()

