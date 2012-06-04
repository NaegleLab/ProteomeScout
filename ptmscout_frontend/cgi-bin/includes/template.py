# None= printHeader(title of page)
# Thus function prints the html code for the top menu of the page and inserts the title argument into the <title> tags
from pathname import *
import sys



import cgi
import MySQLdb
sys.path.append(path+'includes')
sys.path.append(path)
sys.path.append(path+'batch')
sys.path.append(path+'imageScripts')
sys.path.append(path+'dictionaries')
from graph import *
from shapes2 import *
from batchFcns import *
import cgitb; cgitb.enable()
import copy
from sets import Set
from makeDicts import *
import pickle
import re
import cStringIO
import os
import time
import random



def printHeader(title,expid=None,protein_id=None):
	print """<html><head>"""

	if TRACK:
		print """<script type="text/javascript">
		var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
		document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
		</script>
		<script type="text/javascript">
		try {
		var pageTracker = _gat._getTracker("UA-12450206-1");
		pageTracker._trackPageview();
		} catch(err) {}</script>"""
	

	print """<title>%s</title>""" % title
	print """
	<script type="text/javascript" src="%schart/includes/excanvas.js"></script>
	<script type="text/javascript" src="%schart/includes/chart.js"></script>
	<script type="text/javascript" src="%schart/includes/canvaschartpainter.js"></script>
	<script type="text/javascript"  src="%swz_jsgraphics.js"></script>
	<script type="text/javascript"  src="%sgraphs.js"></script>
	<link rel="stylesheet" type="text/css" href="%sptmscout.css">
	<link rel="stylesheet" type="text/css" href="%schart/includes/canvaschart.css" >
	<link rel="stylesheet" type="text/css" href="%stest.css">"""%(jsPath,jsPath,jsPath,jsPath,jsPath,jsPath,jsPath,jsPath)
	

	print """
	<script type="text/javascript">
	function toggle(o) {
	var e = document.getElementById(o);
	e.style.display = e.style.display == 'block' ? 'none' : 'block';
	}
	function toggleGO(d1,d2,d3,form)
	{
	var c = document.getElementById(d1);
	var b = document.getElementById(d2);
	var f = document.getElementById(d3);  
	c.style.display = (c.style.display == 'none') ? 'block' : 'none';
	b.style.display = (b.style.display == 'none') ? 'block' : 'none';
	f.style.display = (f.style.display == 'none') ? 'block' : 'none';
	if (form) {
	if (form.tog.value =="show GO terms")
	    {
	    form.tog.value = "hide GO terms";
	    }
        else
	    {
	    form.tog.value = "show GO terms";
	    }
	}
	}
	"""
	print '</script>'
	print """
	<script type="text/javascript">
	function showmenu(elmnt)
	{
	document.getElementById(elmnt).style.visibility="visible";
	}
	function hidemenu(elmnt)
	{
	document.getElementById(elmnt).style.visibility="hidden";
	}
	</script>
	

	</head>
	<body>"""
	
	print """

	<div id="menu">
	<div style="float:left; text-align:left">
	    <a href="%sindex.py">Experiments</a> |
            <a href='%sbrowse/prosearch.cgi'>Proteins</a>
	  </div>
	  <div style="float:right;text-align:right">
	     <a target='_blank' href="%stitle=Main_Page">Documentation</a> |
	     <a href='%sabout.cgi'>About</a>
          </div>
	</div>

	<div id="page">
	
	<div id="header"><a href="%sindex.py"><img alt="PTMScout logo" src="%slogo.jpg" title="Home"></a></div>

	"""%(urlpath,urlpath,helpPath,urlpath,urlpath,imgPath)
	breadcrumb(expid,title)
	experimentMenu(expid)
	proteinMenu(protein_id)


def loading_indicator(text="Loading"):
	print """<div id="loading">%s</div>""" %(text)

def close_loading_indicator():
	# If threw up a div indicating that data is loading, close it now
	print """<script>document.getElementById("loading").style.display = "none";</script>"""
	
def breadcrumb(expid=None, title=""):
 	print """<div id="breadcrumb">"""
	if expid != None and expid != "None":
		print """<a href="%sexperiment.cgi?expid=%s">Experiment</a> &gt; """%(urlpath,expid)
	print """%s"""%(title)
	print """</div>"""


# Placeholder for displaying a header menu for a selected protein
# call printHeader and provide a protein id as the third argument
def proteinMenu(protein_id=None):
	if (protein_id != None ):
		print '<div id="protein_menu">'
		print '</div>'
	
	
def experimentMenu(expid=None):
	if expid !=None and expid!="None":
		db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
		c=db.cursor()
		menu(c,expid)

def printFooter():
	close_loading_indicator()
	if not TRACK:
		print """<div id="footer">Questions or comments? <br />Email %s
		<div id="footerLicense">
		<p class="box">
		<a rel="license" href="http://creativecommons.org/licenses/by/3.0/"><img alt="Creative Commons License" style="border: medium none ;" src="http://i.creativecommons.org/l/by/3.0/us/88x31.png"/><br /></a>Except where otherwise noted, content on this site is <br /> licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/3.0/">Creative Commons Attribution 3.0 License</a><br>
		<a rel="terms" href=%sterms.cgi>Terms of Use</a>
		</p>
		</div>
		</div>
		</div>
		</body>
		</html>"""%(sysEmail,urlPath)
	else:
		print """<div id="footer">Questions or comments? <br>Email %s
			<div id="footerLicense">
		<p class="box">
		<a rel="license" href="http://creativecommons.org/licenses/by/3.0/"><img alt="Creative Commons License" style="border: medium none ;" src="http://i.creativecommons.org/l/by/3.0/us/88x31.png"/><br /></a>Except where otherwise noted, content on this site is <br /> licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/3.0/">Creative Commons Attribution 3.0 License</a><br>
		<a rel="terms" href=%sterms.cgi>Terms of Use</a>
		</p>
		</div>
		</div>
		</div>
		</body>
		</html>"""%(sysEmail, urlPath)
		
	

# Displays a text-field that is part of a search.  Reads the form
# and displays the current value of the text box.  Arbitrary HTML attributes
# can be passed in with the third argument.
def textboxForSearch(name,form,attribs=""):
	val = form.getvalue(name, "")
	print """<input type="text" name="%s" value="%s" %s />""" %(name,val, attribs)

# Displays a checkbox that is part of a search.  Reads the
# form and decides if the checkbox should have a value or not.
# Arbitrary HTML attributes can be passed in with the third argument.
def checkboxForSearch(name,form,attribs=""):
	if (form.has_key(name)):
		print """<input type="checkbox" name="%s" checked %s />""" %(name,attribs)
	else:
		print """<input type="checkbox" name="%s" %s />""" %(name,attribs)

def radioForSearch(name,val, form,attribs=""):
	fVal = form.getvalue(name, "")
	print """<input type="radio" name ="%s" value="%s" %s""" %(name,val,attribs)	
	if ( fVal == val ):
		print " checked"
	print "</input>"

def selectBoxForSearch(name, values, form,attribs=""):
	fVal = form.getvalue(name,"")
	print """<select name="%s"  %s>"""%(name,attribs)
	if fVal in values:
		values.remove(fVal)
		if fVal !=None:
			values= [fVal]+values
	for item in values:
		if item.strip()!="":
			value = item.strip().replace(" ", "%20")
			print """<option value = '%s'>%s</option>"""%(item,item)
	print "</select>"
	

	
def getExpGroupings(c):
	query = """select id, experiment_id from experiment"""
	c.execute(query)
	x=c.fetchall()
	expDict = {}
	##initialize
	for item in x: expDict[item[1]]=[]
	## dictionary[key=parent]: val=[list of children]
	for item in x:
		if str(item[1])==str(0):
			if item[0] not in expDict.keys(): expDict[item[0]]=[]
		else:
			expDict[item[1]].append(item[0])
	return expDict
		
def makeExpTable(c,form,exp=None,action=urlpath+"index.py"):
    printExpTableJS()
    groupings = getExpGroupings(c)
    
    loading_indicator("Loading Experiments...")

    has_search = False
    if form.has_key("modY") or form.has_key("modT") or form.has_key("modS") or form.has_key("modK") or form.has_key("publishedid") or form.has_key("search_experimentid") or form.has_key("descid") or form.has_key("authorid"):
	    has_search = True

    display_filter = "none"
    display_button = "block"
    if (has_search == True):
	    display_filter = "block"
	    display_button = "none"


    print """
    <div style="height:2em;margin-top:2em;">
    <div style="display:%s;float:left;text-align:left;text-ident:0.2em;">
           <button onClick="document.getElementById('searchExp').style.display='block';this.style.display='none';return true">Filter Experiments</button>
        </div>
	<div style="float-right;text-align:right;margin-right:0.2em;">""" %(display_button)
    print """<a href="%sload/newload.cgi">Load a dataset</a>""" %(urlpath)
    print """
        </div>
      </div>"""
    print """<div style="display:%s;margin-top:1em;margin-bottom:1em;" id="searchExp"><fieldset><legend>Filter Experiments</legend>""" %(display_filter)
    print """
    <form action='"""+action+"""' method="POST">"""
    print """<input type=hidden name='exp_id' value='%s'></td></tr>"""%(form.getvalue("id",exp))
    print """
      <table id="filterExpTable" width="100%" cellspacing=0 cellpadding=0>
      <tr valign="top">
	   <th>Published</th>
	   <th>Experiment</th>
	   <th>Description</th>
	   <th>Author</th>
	   <th>Journal Info</th>
	   <th>Modification</th>
	 </tr>
	 <tr>"""

    print """<td align="center">"""

    checkboxForSearch("publishedid", form)
    
    print "</td><td>"
    textboxForSearch("search_experimentid", form)

    print "</td><td>"
    textboxForSearch("descid", form)
    print "</td><td>"
    textboxForSearch("authorid", form)
    print "<td>"
    checkboxForSearch("modY", form)
    print "Y"
    checkboxForSearch("modT", form)
    print "T<br/>"
    checkboxForSearch("modS", form)
    print "S"
    checkboxForSearch("modK", form)
    print "K"
    print "</td></tr>"
    print """<tr><td colspan=5><input type=hidden name=filter value='true'><input type=submit value='Filter Experiments'>"""
    if (has_search == True):
	    print """<a href="%sindex.py"><button>Clear Search</button></a>""" %(urlpath)
    print "</table></form></fieldset></div>"
    publishedOnly = form.has_key("publishedid")
    Y = form.has_key("modY")
    T = form.has_key("modT")
    S = form.has_key("modS")
    K = form.has_key("modK")
    
    if not True in [Y,T,S,K]:
	    modSearch = False
    else:
	    modSearch =[]
	    if Y: modSearch.append("Y")
	    if T: modSearch.append("T")
	    if S: modSearch.append("S")
	    if K: modSearch.append("K")
    authSearch = form.getvalue("authorid",False)
    descSearch = form.getvalue("descid",False)
    nameSearch = form.getvalue("search_experimentid",False)
    if not(publishedOnly):
        query = """select count(*),experiment.id,name, description, author,primaryModification from experiment join MS on MS.experiment_id=experiment.id group by experiment.id"""
    else:
        query = """select count(*),experiment.id,name, description, author,primaryModification from experiment join MS on MS.experiment_id=experiment.id where published=1 group by experiment.id """
    
    c.execute(query)
    x=list(c.fetchall())
    x.sort()
    x.reverse()
    x = [item[1:] for item in x]
    ### deal with parents
    ## make one row in table for each "parent", then expand to see children

    ###
    if authSearch:
        newx=[]
        for item in x:
		if re.search(authSearch.strip().upper(),item[3].upper()):
			newx.append(item)
        x = newx
    if descSearch:
        newx=[]
        for item in x:
            if re.search(descSearch.strip().upper(),item[2].upper()):
                newx.append(item)
        x = newx
    if nameSearch:
        newx=[]
        for item in x:
            if re.search(nameSearch.strip().upper(),item[1].upper()):
                newx.append(item)
        x = newx
    if modSearch:
	    newx = []
	    for item in x:
		    for res in modSearch:
			    if res in item[4]:
				    newx.append(item)
	    x = newx
    xDict = {}
    for item in x: xDict[item[0]]=item
    string = ''
    string+="""<form action='batchCompare.cgi' method='POST'>"""
    blue = False
    if exp==None:
	    string+="""<table width = "100%"><tr  bgcolor='blue'><th><font color='white'>Experiment</font></th><th><font color='white'>Description</font></th><th><font color='white'>Author</font></th><th><font color='white'>Journal Info</font></th><th><font color='white'>Primary Modifications</font></th></tr>"""
    else:
	    string +="""<input type=hidden value='%s' name='expid'><input type=submit value = 'Compare To Selected'><br><table width = "100%%"><tr  bgcolor='blue'><th></th><th><font color='white'>Experiment</font></th><th><font color='white'>Description</font></th><th><font color='white'>Author</font></th><th><font color='white'>Journal Info</font></th><th><font color='white'>Primary Modification</font></th></tr>"""%(exp)

    ### order the keys
    query = """select experiment_id,count(*) as NUM from MS group by experiment_id order by NUM desc"""
    c.execute(query)
    x=c.fetchall()
    x=[item[0] for item in x]
    orderedGroupings = [item for item in x if item in groupings.keys()]
    for key in orderedGroupings:
        item = xDict.get(key,None)
	if item !=None:
		id, name, desc, author,mod = item
	else: id = None
	if str(id) !=str(exp) and item!=None:
		if len(groupings.get(key,[]))>0:
			buttonText = """<br><a href="javascript:toggleExp('div%s','row%s')"><object id='div%s'>Expand &darr;</object></a><br>"""%(id,id,id)
		else: buttonText=""
		if exp == None:
			string+="""<tr valign='top' onmouseover="this.style.backgroundColor='#cc99ff';" onmouseout="this.style.backgroundColor='%s';" bgcolor="%s"><td><font color = "%s" size="2">"""%(getColor(blue),getColor(blue),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+name+'</a><br>'+buttonText+'</font></td><td><font color = "%s" size="2">'%getFont(blue)+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+desc+'</a></font></td><td><font color="%s" size="2">'%getFont(blue)+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+author+'</a></font></td>'+'<td>%s</td><td><font color="%s" size="2">'%(getCitationString(id,c),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+mod+'</a></font></td>'+'</tr>'
			if len(groupings.get(key,[]))>0:
				counter=0
				for thing in groupings[key]:## make expandable into children experiments
					thing = xDict.get(thing,None)
					if thing!=None: id,name,dex,author,mod = thing
					else: id = None
					if str(id)!=str(exp) and thing!=None:
						string+="""<tr id="row%s_%s" style="display:none" onmouseover="this.style.backgroundColor='#cc99ff';" onmouseout="this.style.backgroundColor='%s';" bgcolor="%s"><td><font color = "%s" size="2">"""%(key,counter,getColor(blue),getColor(blue),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+name+'</a></font></td><td><font color = "%s" size="2">'%getFont(blue)+'<a href="%s'%urlpath+'experiment.cgi?expid=%s">'%(id)+desc+'</a></font></td><td><font color="%s" size="2">'%getFont(blue)+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+author+'<td>%s</td><td><font color="%s" size="2">'%(getCitationString(id,c),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+mod+'</a></font></td>'+'</a></font></td></tr>'
					counter+=1
			
		else:
			string += """<tr valign='top' onmouseover="this.style.backgroundColor='#cc99ff';" onmouseout="this.style.backgroundColor='%s';" bgcolor="%s"><td><input name="include%s"type=checkbox></td><td><font color = "%s" size="2">"""%(getColor(blue),getColor(blue),str(id),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+name+'</a><br>'+buttonText+'</font></td><td><font color = "%s" size="2">'%getFont(blue)+'<a href=%sexperiment.cgi?expid=%s">'%(urlpath,id)+desc+'</a></font></td><td><font color="%s" size="2">'%getFont(blue)+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+author+'<td>%s</td><td><font color="%s" size="2">'%(getCitationString(id,c),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+mod+'</a></font></td>'+'</a></font></td></tr>'
			if len(groupings.get(key,[]))>0:## make expandable into children experiments
				counter=0
				for thing in groupings[key]:## make expandable into children experiments
					
					thing = xDict.get(thing,None)
					if thing!=None: id,name,dex,author,mod = thing
					else: id = None
					if str(id)!=str(exp) and thing!=None:
						string += """<tr id='row%s_%s' style="display:none;" onmouseover="this.style.backgroundColor='#cc99ff';" onmouseout="this.style.backgroundColor='%s';" bgcolor="%s"><td><input name="include%s"type=checkbox></td><td><font color = "%s" size="2">"""%(key,counter,getColor(blue),getColor(blue),str(id),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+name+'</a></font></td><td><font color = "%s" size="2">'%getFont(blue)+'<a href="experiment.cgi?expid=%s">'%id+desc+'</a></font></td><td><font color="%s" size="2">'%getFont(blue)+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+author+'<td>%s</td><td><font color="%s" size="2">'%(getCitationString(id,c),getFont(blue))+'<a href="%sexperiment.cgi?expid=%s">'%(urlpath,id)+mod+'</a></font></td>'+'</a></font></td></tr>'
					counter+=1
		blue = not(blue)
    
    string+='</table></form>'
    string+='<BR><BR>'
    return string

def getCitationString(expid,c):
	query = """select journal,pub_date,volume,pages from experiment where id = %s"""%expid
	c.execute(query)
	x=c.fetchall()
	try:
		journal,pub_date,volume,pages = x[0]
		string = ""
		if journal != "" and journal != None:
			string += journal + ". "
		if pub_date != "" and pub_date != None:
			string += pub_date + ". "
		if volume != "" and volume != None:
			string += "Vol " + str(volume) + ". "
		if pages != "" and pages != None:
			string += pages + "."
	except ValueError: return ""
	except IndexError: return ""
	return string
def getColor(blue):
    if blue: return "#9999cc"
    else: return "white"
def getFont(blue):
    if blue: return "white"
    else: return "black"

def printExpTableJS():
	print"""
	<script type='text/javascript'>
	function toggleExp(pDiv,rowName)
	{
	var counter = 0;
	stillGoing = true;
	
	while(stillGoing==true)
	{
	   if (document.getElementById(rowName+'_'+counter)!=null)
	   { var e = document.getElementById(rowName+'_'+counter);
	     e.style.display = e.style.display == 'table-row' ? 'none' : 'table-row';
	     counter = counter+1;
	     
	     }
	     else {stillGoing = false;}
	}
	
	var f = document.getElementById(pDiv);
	var text = f.innerHTML;
	
	//alert(escape(text).match("u2193"));
	text=escape(text);
	if (text.match("Expand")==null){text=text.replace("Close","Expand");text=text.replace("u2191","u2193");}
	else {text = text.replace("Expand","Close");text=text.replace("u2193","u2191");}
	f.innerHTML=unescape(text);
	
	}
	</script>

	"""


def printExpHeader(expid,c):
    if expid==None: return
    query = """select name, author, date, url, pmid,journal,pub_date,volume,pages,description from experiment where id=%s"""%(expid)
    c.execute(query)
    print """<div id="experimentInfo">"""
    print """<fieldset><legend>Experiment</legend>"""
    try:
	    x=c.fetchall()[0]
	    if len(x)==0: return
	    name, author, date, url, pmid, journal,pub_date,volume,pages,description = x
    except IndexError: return
    if url == "NA": url = None
    if url !=None:
            url = url.replace('&','&amp;')
	    if "http://" not in url: url = "http://" + url
            print """<a href="%s" target="_blank">""" % url

    elif pmid !=None and pmid != 0:
            print """<a href="ww.ncbi.nlm.nih.gov/pubmed/%s" target="_blank">""" % pmid

    if name != "":
	    print name
    if url !=None or pmid !=None:
            print '</a><br>'
    else: print '<br>'
    if author !="":
	    print author + ". "
    if journal != "":
	    print "<b>"+journal + "</b>. "
    if date != "":
	    print pub_date + ". "
    if volume !=None:
	    print "Vol " +str(volume)+". "
    if pages != "":
	    print pages + "."
    
    print "<BR>"
    
    print '<br>%s <br>'% description
    print '<br /> Loaded: %s <br>'% date
    print "</fieldset>"
    print "</div>"



def printPie(fracs,labels,title,type):
	# separate items by colons
	if len(fracs)>10 and len(labels)>10:
		fracs = fracs[0:10]
		labels = labels[0:10]
	fracLabels = str(fracs[0])
	if len(fracs)>1:
		for item in fracs[1:]:
			fracLabels+=':'+str(item)
			
	labelsLabels = str(labels[0])
	if len(labelsLabels)>0:
		for item in labels[1:]:
			labelsLabels += ':'+str(item)
	src = graphDir+"pieGraph.py?type=%s&amp;fractions=%s&amp;labels=%s&amp;title=%s" % (type,fracLabels,labelsLabels,title)
	print """<br><a target="_blank" href="%stitle=Experiment_Summary#Pie_Graphs"><img border="0" alt="help" width="20" src="%shelp.jpg"></a> <a href="javascript:void(window.open('%s','PTMScout Pie Graph','toolbar=no,location=no,directories=no,width=740,height=620,resizable=yes,scrollbars=yes'))">Pie chart of terms</a>"""%(helpPath,imgPath,src)

def helper(title):
	return """<a target="_blank" class="help" href="%stitle=%s"><img class="help" border="0" alt="help" width="20" src="%shelp.jpg"></a>"""%(helpPath,title,imgPath)



strinDict = {'Low':0, 'Medium':1,'High':2}
def printSetForm(table, probeid,multiple=False):
    print '<script type="text/javascript">'
    print 'document.form.table.value = "%s";'%table
    if not(multiple):
        print 'document.form.probeid.value = "%s";'%probeid;
    else:
        for i in range(len(probeid)):
            print 'document.form.elements[%s].value= "%s";'%(str(i+1),probeid[i])
    print '</script>'

def printSetOptions(stringency):
    print '<script type = "text/javascript">'
    print 'document.search.stringency[%s].checked=true;' % strinDict[stringency]
    print '</script>'

def printGOCopyright(c):
	pass

def menu(c,exp_id):
	query = """select id from experiment where export = 0"""
	c.execute(query)
	x=c.fetchall()
	exps = [str(item[0]) for item in x]
	query = """select * from experiment where id = %s"""%exp_id
	c.execute(query)
	x=c.fetchall()
	print "<div id=experimentMenu>"
	print '<a href="%ssummary/summary.cgi?exp_id=%s">Experiment Summary</a> '%(urlpath,str(exp_id))
	if len(x)>0:
		query = """select ambiguity from experiment where id = %s"""%exp_id
		c.execute(query)
		x=c.fetchall()
		try: ambiguity = str(x[0][0])
		except IndexError: ambiguity = "0"
		if str(exp_id) not in exps:
			print '| <a href="%sbrowse/sample.cgi?exp_id=%s">Browse Dataset</a>'%(urlpath,str(exp_id))
			if str(exp_id) not in exps:
				print '| <a href="%ssubset/select.cgi?exp_id=%s">Evaluate Subsets</a>'%(urlpath,str(exp_id))
				print '| <a href="%sbatch/batch.py?exp_id=%s">Compare Datasets</a>'% (urlpath,str(exp_id))
			
				if ambiguity == "1":
					print "| <a href='%sambiguity/ambiguity.cgi?expid=%s'>Report Ambiguity</a>"%(urlpath,exp_id)
	print """</div>"""

def printGODocs(c):
	print "<a target='_blank' href='%stitle=References_to_Data_Sources#Gene_Ontology'>Gene Ontologies (GO)</a> %s"%(helpPath,helper("References_to_Data_Sources#Gene_Ontology"))
	# get version number and date
	date=''
	ver=''
	query = """select version from GO"""
	c.execute(query)
	x=c.fetchall()
	if len(x)>0: ver = x[0][0]
	query = """select version from protein_GO"""
	c.execute(query)
	x=c.fetchall()
	if len(x)>0: date = x[0][0]
	print "<BR><font size=3>Version %s %s</font><BR>"%(ver,date)
	


def getAllSpecies(c):
	query = """select species from protein group by species"""
	c.execute(query)
	x = c.fetchall()
	return [item[0] for item in x]



def allSameExperiment(msids,c):
	exps = []
	for item in msids:
		query = """select experiment_id from MS where id = %s"""%item
		c.execute(query)
		x=c.fetchall()
		try:
			if int(x[0][0]) not in exps:
				exps.append(int(x[0][0]))
		except IndexError:
			pass
	if len(exps) > 1:
		return False
	else:
		return True
		
			   
