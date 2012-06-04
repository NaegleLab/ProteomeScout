from pathname import *
import sys
sys.path.append(path+'includes')
sys.path.append(path)
from template import *
from peptideFunctions import *
from tableFunctions import *
import copy
def getMSids_batch(peptides,c,expid):
    """gets the MS ids matching a group of phosphopep_ids from a specific experiment_id"""
    msids = []
    for peptide in peptides:
        query = """select MS_id from MS_phosphopep join MS on MS_phosphopep.MS_id = MS.id join protein on protein.id=MS.protein_id where phosphopep_id = %s and experiment_id=%s order by acc_gene"""%(peptide,expid)
        c.execute(query)
        x=c.fetchall()
        x=[item[0] for item in x]
        msids = msids + x
    return msids
def makeTable(difference,c,expid,formdata,letter,form):
    """ makes table of novel peptides"""
    newms = alphabetizeMS(getMSids_batch(difference,c,expid),c)
    print '<table width="60%%"><tr><td>'
    print '<BR>'
    if len(newms)>0:
        #exp_id, c,search,stringency,form,cols,url,letter=None,msids=[],width='100',syn=True,meta_queries=[]
        #makeWholeTable(newms, c,'Low',expid,letter=letter,user=user)
        makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","pep_aligned"],"",msids=newms,width='70',newtab=True)
        print '<BR><BR><BR>'
    print '</td></tr></table>'




def makeExpTable(c):
    query = """select id,name, description, author from experiment"""
    c.execute(query)
    x=c.fetchall()
    string = ''
    string+='<table width = "100%">'
    blue = False
    string+="""<tr  bgcolor='blue'><th><font color='white'>Include</font></th><th><font color='white'>Experiment</font></th><th><font color='white'>Description</font></th><th><font color='white'>Author</font></th></tr>"""
    for item in x:
        id, name, desc, author = item
        string+="""<tr onmouseover="this.style.backgroundColor='#cc99ff';" onmouseout="this.style.backgroundColor='%s';" bgcolor="%s"><td><input name="include%s"type=checkbox></td><td><font color = "%s" size="2">"""%(getColor(blue),getColor(blue),str(id),getFont(blue))+'<a href="experiment.cgi?expid=%s">'%id+name+'</a></font></td><td><font color = "%s" size="2">'%getFont(blue)+'<a href="experiment.cgi?expid=%s">'%id+desc+'</a></font></td><td><font color="%s" size="2">'%getFont(blue)+'<a href="experiment.cgi?expid=%s">'%id+author+'</a></font></td></tr>'
        blue = not(blue)
    
    string+='</table>'
    string+='<BR><BR>'
    return string
def getColor(blue):
    if blue: return "#9999cc"
    else: return "white"
def getFont(blue):
    if blue: return "white"
    else: return "black"

## get all possible phosphopeps something could be
def getPossiblePeps(sequence,c):
    sequence = sequence.upper()
    proteins = getProteins(sequence,c)
    possible = []
    for protein in proteins:
        query = """select id,pep_aligned from phosphopep where protein_id = %s""" %protein[1]
        c.execute(query)
        x = c.fetchall()
        for item in x:
            id,seq = item
            seq = seq.upper()
            if seq in sequence: possible.append(id)
    return possible
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
def alphabetizeMS(msids,c):
    """return the msids ordered by acc gene"""
    ms = []
    for item in msids:
        query = """select acc_gene from protein join MS on MS.protein_id=protein.id where MS.id = %s"""%item
        c.execute(query)
        x=c.fetchall()
        if len(x)>0:
            ms.append((x[0][0],item))
        else: ms.append(('',item))
    ms.sort()
    return [item[1] for item in ms]
