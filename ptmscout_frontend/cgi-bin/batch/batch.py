#!/usr/local/bin/python
from pathname import *
import cgi
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import sys

sys.path.append(path+'includes')
from template import *
def main():
    db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    c=db.cursor()
    form = cgi.FieldStorage()
    try:
        exp_id = int(form.getvalue('exp_id',None))
    except TypeError:
        printHeader("Compare Datasets", None)
        print "<BR><BR>No experiment selected<BR><BR>"
        printFooter()
        sys.exit(0)
    letter = form.getvalue("letter","ALL")

    printHeader("Compare Datasets", exp_id)
    printJS()
    
    expTest = int(form.getvalue('expTest',0))
    if expTest == 0:
        expTests = []
    else:
        query = """select name from experiment where id = %s"""%str(expTest)
           ## get nameExp
        c.execute(query)
        x=c.fetchall()
        nameExp = x[0][0]
        expTests = [(expTest,nameExp)]
    printExpHeader(exp_id,c)
   
    print '<a name="home"></a>'

    loading_indicator();
    print """<BR><BR><form method="POST" action="novel.cgi"><input type=hidden value=%s name="expid"><table><tr><td>"""%exp_id
    print """<input type=submit name = "novel" value = "View novel sites"></td><td>List all sites found only in this experiment %s</td></tr><tr><td></td><td><input type=checkbox name="published">Include only published datasets<br><br></td></tr>"""%(helper("Compare_With_Standard_Datasets#Novel_Sites"))
    print """<tr><td><input type=button name = "compare" value = "Compare Experiments" onclick="tog('expTable');"></td><td>Compare sites across multiple experiments %s </td></tr>"""% helper("Compare_With_Standard_Datasets#Across_Experiments")
    print "</table></form>"
    if form.has_key("filter"):
        print """<BR><BR><div style="width:70%;display:block" id="expTable">"""
    else:
        print """<BR><BR><div style="width:70%;display:none" id="expTable">"""
    print makeExpTable(c,form,exp=exp_id,action="batch.py")
   
    print "</div>"
    
    printFooter()
   
    
def getMSids(peptides,c,expid):
    msids = []
    for peptide in peptides:
        query = """select MS_id from MS_phosphopep join MS on MS_phosphopep.MS_id = MS.id where phosphopep_id = %s and experiment_id=%s"""%(peptide,expid)
        c.execute(query)
        x=c.fetchall()
        x=[item[0] for item in x]
        msids = msids + x
    return msids


def printJS():
    print """
    <script type="text/javascript">
	
	function tog(o)
	{
	    var e = document.getElementById(o);
	    e.style.display = e.style.display == 'block' ? 'none' : 'block';   
	    
	}
	</script>
	"""

    
main()
