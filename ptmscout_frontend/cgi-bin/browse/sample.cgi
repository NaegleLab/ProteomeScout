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
        search = form.getvalue("acc_search","no search")
	search_display = form.getvalue("acc_search", "")
	letter=form.getvalue("letter","A")
	stringency = form.getvalue("stringency","Low")
	# if I knew pythong better, this would be so bad
	low_checked = ""
	med_checked = ""
	high_checked = ""
	if stringency == "High":
		high_checked = "checked"
	elif stringency == "Med":
		med_checked = "checked"
	else:
		low_checked = "checked"
	 
	exp_id= form.getvalue("exp_id",None)
	printHeader("Browse Dataset",expid=exp_id)
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
	c=db.cursor()
	query = """select id from experiment where export = 0"""
	c.execute(query)
	x=c.fetchall()
	exps = [str(item[0]) for item in x]
	if exp_id in exps:
		print "<BR><BR>Illegal experiment id<BR><BR>"
		exp_id=None
	else:
		printExpHeader(exp_id,c)
	
		print """
		<fieldset class="filter"><legend>Dataset Filter</legend>
	        <form name ="search" method="POST" action ="sample.cgi?exp_id=%s">		
		<input type = "hidden" name = "letter" value = "%s" />
		<input type="hidden" name="exp_id" value = "%s" />
	        Protein Name: <input type="text" name="acc_search" size="50" value="%s" />
		 <br />Scansite Stringency:
		 <div style="margin-left:2em;margin-right:2em;">
  		 <input type="radio" name=stringency value = "Low" %s /> Low <br/>
		 <input type="radio" name = stringency value = "Medium" %s /> Medium <br />
	         <input type="radio" name = stringency value = "High" %s /> High
		 </div>
		 <div style="margin-top:1em;text-align:left"><input type="submit" value="Filter" /></div>
	</fieldset></form>""" % (exp_id,letter,exp_id,search_display, low_checked, med_checked, high_checked )
	        printSetOptions(stringency)
		if exp_id!=None:
			makeWholeExpTable(exp_id,c,search,stringency,form,["MSid","sequence","protein","species","site","pep_aligned","predictions"],"sample.cgi",letter=letter,width='75%')
		else:
			print "<h3>No experiment selected</h3>"
	print '<br><br>'
	db.close()
	printFooter()








	
main()

