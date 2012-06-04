#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import sys

sys.path.append(path+'includes')
from template import *
printHeader("PTMScout Error")
#print """ <div id="content">"""
print "<html>"
print "<b> We're sorry, an error has occured. <BR> Please contact ptmscout_admin@mit.edu with information regarding the experiment you were working on (if any) and the operations you were performing.  <BR>An estimate of the day and time will also be helpful in debugging."
print "</html>"
#print "</div>"
#db.close()
printFooter()


