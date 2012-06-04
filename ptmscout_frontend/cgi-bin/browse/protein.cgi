#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import sys

sys.path.append(path+'includes')
from template import *
import checks
from proteinFcns import *

# Prints relative links to the various sections on this long page
# XXX BOGON ALERT - this is basically duplicate in ambiguity/compare.cgi
def startSection(title, anchor, help=""):
	print """<h2><a name="%s">%s %s</a></h2>""" %(anchor, title, helper(help))
	printNavigation(anchor)
	print '<div class="section">'

def endSection():
	print '</div>'

def printNavSection(anchor, title, display_location,last=False):
	if ( anchor == display_location ):
		print title
	else:
		print """<a href="#%s">%s</a>""" %(anchor, title)
	if (last != True ):
		print '|'
		

def printNavigation(display_location=""):
	print '<div class="protein_navigation">'
	printNavSection('details'	, 'Protein Details'		, display_location)
	printNavSection('gene_ontology'	, 'Gene Ontologies'		, display_location)
	printNavSection('graphs'	, 'Expression Data'		, display_location)
	printNavSection('domains'	, 'Domain Structure'		, display_location)
	printNavSection('mods'		, 'Modification Sites'	, display_location)
	printNavSection('data'		, 'Data'			, display_location, True)
	print "%s" %(helper("Protein_Page"))
	print '</div>'
# END BOGON ALERT

# The main rendering function for this page
def main():
	### get form data
	form = cgi.FieldStorage() 
	pid = form.getvalue("protein_id","3810")
	exp_id=form.getvalue("exp_id",None)
	printHeader("Protein Details",expid=exp_id,protein_id=pid)
	checks.printJS()
	try:
		referrer = cgi.os.environ["HTTP_REFERER"]
	except KeyError:
		referrer = "None"
	if "ambiguity" in referrer or "subset" in referrer or "batch" in referrer: referrer = None

	### connect to database
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
	c=db.cursor()
	### end connect to database
	if exp_id!=None and exp_id!="None":
		printExpHeader(exp_id,c)


	prodict = makeProteinInfoDict(pid,c)
	
	print '<div id="content">'
	makeProteinDetails(prodict)
	makeGeneOntologies(prodict,c)
	makeGraphs(form,prodict,pid,c,exp_id)
	makeDomainStructure(pid,c,exp_id)
	makePhosphorylationSites(pid,c)
	makeData(pid,exp_id,c)
	print '</div>'
	
       	db.close()
	printFooter()


## get which phoshopeps to graph, return a list of tuples [(MS_id, label)] (in the form Y234 etc.)
def getMStoGraph(pid,exp_id,c):
	c.execute("""select id from MS where protein_id=%s and experiment_id=%s""",(pid,exp_id))
	x=c.fetchall()
	msids = [item[0] for item in x]
	names =[]
	for ms in msids:
		c.execute("""select site_type,site_pos from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id = %s""",(ms,))
		x=c.fetchall()
		names.append(x[0][0]+str(x[0][1]))
	return [(msids[i], names[i]) for i in range(len(msids))]

	

def makeData(pid,exp_id,c):
	chart = 'bar'
	errors = 'true'
	numSites = len(getMStoGraph(pid,exp_id,c))
	       	
	# first see if there is actually any data
	c.execute("""select * from data join MS on data.MS_id=MS.id where experiment_id=%s and protein_id=%s""",(exp_id,pid))
	x=c.fetchall() #x will be checked later
	
	startSection("Data", 'data', "Protein_Page#Data_Plots")

	#### start printing data
	## make javascript function to toggle errorbars
	print '<script type="text/javascript">'
	print 'function redraw() {'
	print 'var d = document.getElementById("graph")'
	print 'var transform = document.graphform.transform.value;'
	print 'var check = document.graphform.errorbars.checked'
	print ' if (check == false)' # if errors are true
	print """{d.innerHTML = "<img  alt='graph unavailable' src='%smakeGraph.py?exp_id=%s&pid=%s&errors=%s&amp;chart=%s&amp;transform="+transform+"'>"} """ % (graphDir,exp_id,pid,errors,chart)
	
	print 'else'
	print """ {d.innerHTML = "<img  alt='graph unavailable' src='%smakeGraph.py?exp_id=%s&pid=%s&errors=%s&transform="+transform+">"}  """ % (graphDir,exp_id,pid,'true')
	print '}' # end redraw function
	print '</script>'
	
	if len(x)>0:   	
		print """<form name = graphform method = GET action = ""><input type = checkbox name = 'errorbars' onChange="redraw()" checked>Include Errobars<br>"""
		if numSites>1:
			print '<p>Group by:&nbsp;&nbsp; <select name = group onChange="redraw()"><option value = run>Run</option><option value = site>Site</option></select>'
		print '<p>Select Data Transform: &nbsp;&nbsp;<select name = transform onChange = "redraw()"><option value = None>None</option><option value = normToMax>Normalize to Maximum</option>'
		print '</select>'
		print '</form>'
		print '<div id = "graph">'
		print '<img alt="graph" src="%smakeGraph.py?exp_id=%s&amp;pid=%s&amp;errors=%s&amp;chart=%s">' % (graphDir,exp_id,pid,errors,chart)
		
		print '<a name = "drawGraph"></a>'
		print '</div>'

	print "<BR><a target = '_blank' href='site.cgi?pid=%s'>View all data for this protein.</a><BR><BR><BR>"%pid
	endSection()

def getPeptide(c,item):
	query = """select pep_aligned from phosphopep where id = %s"""%item
	c.execute(query)
	try:
		pep = c.fetchall()[0][0]
	except IndexError: pep = ''
	return pep
	
	
# Draws the phosphorylation section of the page
def makePhosphorylationSites(pid,c):
	startSection("Modification Sites", "mods", "Protein_Page#Site_Experiment_Map")

	# Create the dictionary of data we'll draw the tables
	pDict={}
	c.execute("""select id from phosphopep where protein_id=%s order by site_pos""",(pid,))
	x=c.fetchall()
	x=[thing[0] for thing in x]
	for item in x:
		info={}
		c.execute("""select site_type,site_pos from phosphopep where id=%s""",(item,))
		id=c.fetchall()[0]
		name = id[0]+str(id[1])
		info["name"]=name
		c.execute("""select experiment_id from MS join MS_phosphopep on MS_phosphopep.MS_id=MS.id where phosphopep_id=%s""",(item,))
		exp=c.fetchall()
		exp=[e[0] for e in exp]
		#uniquify
		exp = reduce(lambda l, x: x not in l and l.append(x) or l, exp, [])
		# exp_info is a list of tuples of name, author of experiment
		exp_info=[]
		for e in exp:
			c.execute("""select name, author, URL, PMID, id from experiment where id = %s""",(e,))
			y=c.fetchall()[0]
			exp_info.append((y[0],y[1], y[2],y[3],y[4]))
		info["exp"]=exp_info
		# pDict is a dictionary with key=experiment id, value is tuple of name, author, url, pmid
		pDict[item]=info


	print '<table id="peptide_table" cellspacing=0 cellpadding=0>'
	
	print '<tr><th>Site</th><th>Peptide</th><th>Experiment</th></tr>'

	count = 0
	for item in x:
		count = count + 1
		style = "odd"
		if ((count % 2) == 0):
			style = "even"
		print """<tr valign="top" class="%s">""" %(style)

		# Site
		name = pDict[item]["name"]
		print """<td><a name="%s">%s</a></td>""" %(name,name)

		# Peptite
		print '<td>'
		print getPeptide(c,item)
		print '</td>'
		
		# table with experiment data
		print '<td>'
		print '<ul class="no_bullet">'
		for i in range(len(pDict[item]["exp"])):
 			name,author, url, pmid,id=pDict[item]["exp"][i][0],pDict[item]["exp"][i][1],pDict[item]["exp"][i][2],pDict[item]["exp"][i][3],pDict[item]["exp"][i][4]
			
 			print """<li><a href="%sexperiment.cgi?expid=%d">%s-%s</a></li>"""%(urlpath,id,name,author)
# 			if url !=None and url !='UNPUBLISHED':
# 				url = url.replace('&','&amp;')
# 				url = url.replace('-','&ndash;')
# 				newname=''
# 				for i in range(len(name)):
					
# 					if name[i].isalnum() or name[i] ==' ':
# 						newname += (name[i])
# 				name = newname
					
				
# 				print '<li><a href="' + url+'" target="_blank">'+name+' - ' + author + '</a></li>'
# 			elif pmid != None:
# 				print '<li><a href="' + 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&amp;cmd=search&amp;term=' + str(pmid) +'"target="_blank">'+name+' - ' + author + '</a></li>'
# 			else:
# 				print '<li>'+name+' - ' + author + '</li>'
		print '</ul>'		
		print '</td></tr>'
		
	print '</table>'
	endSection()

# Draws the domain structure
def makeDomainStructure(pid,c,exp_id):
	startSection("Protein Structure", "domains", "Protein_Page#Protein_Structure")	

	d=makeDict(pid,c)
	if '~~~' in d: del d['~~~']
	p=makepDict(pid,c)
	olist,elist = makePOthersList(pid,c,exp_id)
	
		
	print "<div id='Canvas' style='position:relative;bottom:150px;'></div>"
	
	length = getLength(pid,c)
	p = Protein(pid,length,d,p,elist,olist,c)
	p.drawInJS(exp_id)

	endSection()
	

# Helper for making graphs
def printForm(pid,species,c,exp_id):
	data='tissue'
	probes = getProbeOptions(pid,c)
	print '<div style="background-color:#ccccff;border:1px gray solid;width:100%;padding:3px;">'
		# fill in action
	print '<form name = "form" method=POST action="%sbrowse/protein.cgi?exp_id=%s#graphs">' % (urlpath,str(exp_id))
	print '<input type="hidden" name="protein_id" value="%s">' % pid
	print '<input type="hidden" name="exp_id" value = "%s">' % str(exp_id)
	print '<p align=left>'
	print 'Probeset ID:'
	print '<select name= "probeid">'
	# one of these for each probe - make for loop#
	for probe in probes:
		print '<option value = "%s">%s</option>' % (probe,probe) # fill in with probes
	print '</select>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
	print 'Data Type:'
	print '<select name = "table" onChange="changeChecks(this.form)">'
	if species =="homo sapiens":
		print '<option value="cell">NCI60</option>'
		print '<option value="tissue">Human Expression</option>'
		print '<option value="both">View Both</option>'
		if data =='tissue': cols=EXPRESSION_HUMAN_COLS
		else: cols=EXPRESSION_NCI60_COLS
	else:
		print '<option value = "null">Mouse Expression</option>'
		cols = EXPRESSION_MOUSE_COLS
	print '</select>'
	print '<p align=left>'
	print 'Select columns to display: <input type=button name="buttons" value="show column choices" onClick="addRemoveChecks(this.form)"><br>'
	print '<INPUT TYPE="submit" NAME="button" Value="Display">'
	 

	# begin checkbox div
	print '<div id="checks">'
	print '</div>'
	print '</form>'
	print '</div>'


def makeGraphs(form,prodict,pid,c,exp_id):
	# XXX BOGON ALERT - This seems duplicate in ambiguity/compare.cgi
 	table 	= form.getvalue("table","")
	columns = form.getvalue("columns","")
	probeid	= form.getvalue("probeid","")

	ht,hc,mo = [],[],[]
	if form.has_key("human_tissue_cols"):
		ht = form.getlist("human_tissue_cols")
	if form.has_key("nci60_cols"):
		hc = form.getlist("nci60_cols")
	if form.has_key("mouse"):
		mo = form.getlist("mouse")
	if table=="tissue":
		columns = ht
	elif table=="cell":
		columns = hc
	elif table =="mouse":
		columns = mo
	else:
		columns = []

	if probeid=="":
		x=getProbeOptions(pid,c)
		if len(x)>0: probeid = x[0]
		else: probeid = ''
	if table =="":
		table = 'tissue'
	if columns==[] or 'all' in columns:
		columns = 'all'
	# END BOGON ALERT

	startSection("Expression Data", "graphs", "Expression_Data")
	
	print "<h6>Expression Data comes from the Genomics Institute Of The Novartis Research Institute <a target='_blank' href='http://symatlas.gnf.org/'>SymAtlas Project</a></h6>"
	print "<h6><a href='%sdisclaimer.html' target='_blank'>GNF SymAtlas Copyright</a></h6>"%jsPath
	if probeid!='':
		printForm(pid,prodict["species"],c,exp_id)


	printSetForm(table,probeid)
	if table=='both' and probeid != '':
		columnsBoth= [ht,hc]
		if ht==[]:
			columnsBoth[0]='all'
		if hc==[]:
			columnsBoth[1]='all'
		makeGraph(probeid,prodict["species"],c,table='tissue',columns=columnsBoth[0])
		makeGraph(probeid,prodict["species"],c,table='cell',columns=columnsBoth[1])
	elif probeid!='':
		makeGraph(probeid,prodict["species"],c,table=table,columns=columns)
	else:
		print '<h3>No probes found</h3>'
		
	endSection()



# Displays the detailed information about the protein
def makeProteinDetails(dictionary):
	url_dict = makeUrlDict()	

	startSection("Protein Details - " + dictionary["acc"], "details")
	print '<p><label>Name:</label>'
	print dictionary["name"]
	print '</p>'
	print '<p><label>Species:</label>'
	print dictionary["species"]
	print '</p>'
	print '<p><label>Accessions:</label>'
	print '<ul>'
	for item in dictionary["accs"]:
		print '<li><i>'
		print item[1] + '</i>: '
		if item[1]!='gene_synonym':
			print """<a href='%s""" % url_dict.get(item[1])
			print """%s' target='_blank'>""" % item[2]
		else:
			print str(item[2])
		if item[1]!='gene_synonym':
			print """%s</a>""" % item[2]
		print "</li>"
	
	print '</ul>'
	print '</p>'
	endSection()
	
def makeGeneOntologies(dictionary,c):
	url_dict = makeUrlDict()	
	go_dict = makeGODict()

	startSection("Gene Ontologies - " + dictionary["acc"], "gene_ontology")

	print '<div style="margin-bottom:1em">'
	printGODocs(c)
	print '</div>'
	
	print '<table cellspacing=0 cellpadding=0 id="gene_ontology_breakdown">'
	print '<tr>'
	print '<th>Molecular Functions</th>'
	print '<th>Cellular Components</th>'
	print '<th>Biological Process</th>'
	print '</tr>'
	
	
	print '<tr>'	
	print '<td valign="top">'
	makeGOList('F', dictionary, url_dict, go_dict)	
	print '</td>'

	print '<td valign="top">'
	makeGOList('C', dictionary, url_dict, go_dict)	
	print '</td>'

	print '<td valign="top">'
	makeGOList('P', dictionary, url_dict, go_dict)	
	print '</td>'
	print '</tr>'
	print '</table>'

	endSection()
	

	


				
	
main()
