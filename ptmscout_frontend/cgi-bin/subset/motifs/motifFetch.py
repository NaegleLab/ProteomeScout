#!/usr/bin/python
import random
import cgi
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
from pathname import *
form = cgi.FieldStorage()
num = form.getvalue('num',80209)
if num != None:
    try:
        f = open(motifPath+'htmls/outHTML'+str(num),'r')
        text = f.read()
        print text
        f.close()
    except IOError: print "<H3>Error retrieving results. Results might not be ready yet. Please try again later</h3>"
else:
    print "Error"
