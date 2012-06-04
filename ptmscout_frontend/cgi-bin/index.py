#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import sys

sys.path.append(path+'includes')
from template import *
db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
c=db.cursor()
form = cgi.FieldStorage()
printHeader("Home - Experiments")
print """ <div id="content">"""
print """   <table width = "100%"><tr><td align="justify">PTMScout is a web application for viewing and analyzing data from mass spectrometry experiments and other datasets for post-translational protein modifications. The goal of PTMScout is to provide a useable interface for biologists to easily work with large datasets of post-translational measurements.</td></tr></table>"""
print makeExpTable(c,form)
print "</div>"
db.close()
printFooter()


