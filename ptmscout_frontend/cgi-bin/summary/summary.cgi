#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import sys

sys.path.append(path+'includes')
from template import *
try:
	f = open(path+'dictionaries/GOdict.dict.pkl','rb')
	GODict = pickle.load(f)
	f.close()
except IOError: GODict = {}

def main():
	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
	c=db.cursor()
	form = cgi.FieldStorage()
	exp_id= form.getvalue("exp_id",None)
	if exp_id == None:
		print "Error, null experiment id. </body></html>"
		sys.exit(0)
	printHeader("Experiment Summary",exp_id)
	exp_id=int(exp_id)
	query = """select id from experiment"""
	c.execute(query)
	x=c.fetchall()
	if exp_id not in [item[0] for item in x]:
		print "Error, invalid experiment id.</body></html>"
		sys.exit(0)
	
	query = """select export from experiment where id = %s"""%exp_id
	c.execute(query)
	x=c.fetchall()
	try:
		export = bool(int(x[0][0]))
	except IndexError: export = False
	printExpHeader(exp_id,c)
	loading_indicator("Loading ...")
	try:
		f = open(summaryPath+'exp%s.txt'%exp_id,'w')
	except IOError:
		f = open(summaryPath+'exp%s.txt%s'%(exp_id,random.randint(0,1000)),'w')
	### print number of each thing
	files = os.listdir(path+"dictionaries/savedDicts")
	exps = []
	for item in files:
		try: exps.append(int(item[4:]))
		except ValueError: pass
	
	printSummary(c,exp_id,f)
	print """<a href='summary.txt?filename=%s'>Export summary to text file</a>  %s<BR><BR>"""%(f.name,helper("Experiment_Summary#Export_to_Text_File"))
	printJavaScript()		
	allmsids = getMSids(exp_id,c)
	numMS,numPro,numPeps,MS, Peps,Pro = getParameters(c,getMSids(exp_id,c))
	qs = "numMS=%s&numPro=%s&numPep=%s&select=all&exp_id=%s&dataType=metadata&strin=medium&stringency=medium"%(numMS,numPro,numPeps,exp_id)
	### get exps we have dictionaries for
	
	if len(allmsids)>0:
		
		msids = []
		
		if numPeps<4000 or int(exp_id) in exps:
			if int(exp_id) in exps:
				import pickle
				file= open('%sdictionaries/savedDicts/dict%s'%(path,exp_id),'r')
				dall = pickle.load(file)
				dsub = dall
				file.close()
			else:
				dall = makeSummaryDict(exp_id,c)
				dsub = dall
			print '<br><font size=5><b>GO Data</b> </font> %s<BR><BR>'%(helper("Experiment_Summary#GO_Terms"))
			print '<div id = "GO" style="display:block">'
			printGODocs(c)
			printGODiv(dsub,dall,msids,allmsids,f,qs,expid=exp_id,export=export)
			print '</div>'

			print '<br><font size=5><b>Pfam Site</b></font> %s<BR><BR>'%(helper("Experiment_Summary#Pfam_Sites"))
			print """<input type=button name = PFAMButton value = "Show PFAM" onclick="tog('PFAM','PFAMButton');">"""
			print '<div id = "PFAM" style="display:none">'
			printPfamDiv(dsub,dall,msids,allmsids,f,qs,expid=exp_id,export =export)
			print '</div>'

			print '<br><font size=5><b>Pfam Domains</b></font> %s<BR><BR>'%(helper("Experiment_Summary#Domains"))
			print """<input type=button name = DOMButton value = "Show Domains" onclick="tog('DOM','DOMButton');">"""
			print '<div id = "DOM" style="display:none">'
			printDomDiv(dsub,dall,msids,allmsids,f,qs,expid=exp_id,export= export)
			print '</div>'

			print '<br><font size=5><b>Predictions</b></font> %s<BR><BR>'%(helper("Experiment_Summary#Scansite_Predictions"))
			print """<input type=button name = SCANButton value = "Show Scansite" onclick="tog('SCAN','SCANButton');">"""
			print '<div id = "SCAN" style="display:none">'
			printScanDiv(dsub,dall,msids,allmsids,f,qs,expid=exp_id,export= export)
			print '</div>'
			print '<br><Br><br>'


	else:
		print 'No Data Found<BR><BR>'
	db.close()
	f.close()
	printFooter()
	
def printSummary(c,expid,f):
	numMS, numPro, numPeps = getExperimentParameters(c,expid)
	c.execute("""select name, author, date, url, pmid from experiment where id=%s""",(expid,))
	x=c.fetchall()[0]
	name, author, date, url, pmid = x
	f.write("Name:\t"+name+"\nAuthor:\t"+author+"\nDate:\t"+date+"\n")
	sites = getSiteTypes(expid,c)
	species = getSiteBySpecies(expid,c)
	print '<hr  width = "70%%">'
	print '<table width = 600>'
	print """<tr><td>%s<BR></td></tr>"""%(helper("Experiment_Summary#Experiment_Overview"))
	print '<tr><td>Measured Peptides</td><td>%s</td></tr>' % numMS
	print '<tr><td>Proteins</td><td>%s</td></tr>' % numPro
	print '<tr valign = top><td>Modification Sites</td><td>%s</td></tr>'% numPeps
	print "<tr valign = top><td>"
	print "<table>"
	print "<tr><td>&#8226; By site type:</td></tr>"	
	for item in sites:
		print "<tr><td>%s:</td><td>%s</td></tr>"%(item, sites[item])
	
	print "</table>"
	print "</td><td><td></td></tr><tr><td>"
	print "<table cellspacing = '2' cellpadding = '2'>"
	print "<tr><td><br>&#8226; By species:</td></tr>"
	for item in sortListByElement(species.items(),1):
		print "<tr><td>%s:</td><td>%s</td></tr>"%(item[0], item[1])
		print "<tr><td><table>"
		for s in getTypeBySpecies(item[0],expid,c).items():
			print "<tr><td>%s:</td><td>%s</td></tr>"%(s[0],s[1])
		print "</table></td></tr>"
	print "<tr><td><BR><BR></td><td><br><BR></td></tr>"
	print "<tr><td>Rejected peptides on experiment load</td><td>%s</td></tr>"%getNumRejected(expid,c)
	print "<tr><td>%s</td><td></td></tr>"%('<a href="'+urlPath+'load/errors.cgi?expid=%s" target="_blank">View Error Log</a>'%expid)
	#print "<tr><td><a href='%ssummary/exportDataset.txt?expid=%s'>View entire dataset</a></td><td></td></tr>"%(urlPath,expid)
	print "<tr><td><BR><BR></td><td><br><BR></td></tr>"
	print "</table></td></tr>"
	
	###make weblogo file

	filename = "%sexp%s.fasta"%(fastaPath,expid)
	makeFasta(c,expid,filename)
	print "<tr><td>Peptide sequence profile: %s</td></tr>"%helper("Weblogo")
	print "<tr><td><img  width='400' src='%sweblogo.cgi?filename=%s'></td></tr>"%(graphDir,filename)
	print '</table>'
	
	f.write("Peptides:\t"+str(numMS)+"\nProteins:\t"+str(numPro)+"\nPhosphorylation Sites:\t"+str(numPeps)+"\n")
	print '<hr  width = "70%%">'
	f.write('\nSites by type:\n')
	for item in sites:
		f.write(item+':'+str(sites[item])+'\n')
	f.write('\nSites by species:\n')
	for item in species:
		f.write(item+':'+str(species[item])+'\n\n')

def getSiteBySpecies(expid,c):
	query = """select species, count(*) from phosphopep join MS_phosphopep on phosphopep.id = MS_phosphopep.phosphopep_id join MS on MS_phosphopep.MS_id=MS.id join protein on MS.protein_id=protein.id where experiment_id=%s  group by species"""%expid
	c.execute(query)
	x=c.fetchall()
	return dict(x)
	
def getSiteTypes(expid,c):
	query = """select site_type,count(*) from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id = phosphopep.id join MS on MS.id = MS_phosphopep.MS_id where experiment_id=%s group by site_type"""%expid
	c.execute(query)
	x=c.fetchall()
	newx = []
	for item in x:
		new = ''
		if item[0]=="Y": new = "Tyrosine"
		if item[0]=="S": new = "Serine"
		if item[0] == "T":new = "Threonine"
		if item[0]=="K": new = "Lysine"
		newx.append((new,item[1]))
	return dict(newx)
def getExperimentParameters(c,expid):
    ## get number of MSids
    query = "select * from MS where experiment_id=%s" % expid
    c.execute(query)
    x=c.fetchall()
    numMS = len(x)
    ## get number of proteins
    proteins = sets.Set([item[3] for item in x])
    numPro = len(proteins)
    ## get number of peptides
    query = "select * from MS_phosphopep join MS on MS_phosphopep.MS_id = MS.id where experiment_id=%s"% str(expid)
    c.execute(query)
    x=c.fetchall()
    numPeps = len(x)
    return numMS, numPro,numPeps




	

def printGODiv(dsub,dall,msids,allmsids,f,qs,sort="enrichment",expid=None,export = True):
	dsub = dall
	msids = allmsids
	f.write(">GO: Biological Process\n")
	print """<BR><input type=button value = 'Show GO:Biological Process' name = 'GO_BPDivButton' onclick = 'tog("GO_BPDiv","GO_BPDivButton")'>"""
	print """<div id = 'GO_BPDiv' style="display:none">"""
	print '<center><h3>GO: Biological Process</h3></center>'
	subsetGOTable(dsub['GO_BP'],dall['GO_BP'],len(allmsids),len(msids),f,sort="enrichment",linkType="GO",qs=qs,heading="GO_BP",expid=expid,export = export)
	
	print "<br><br>"
	x=getGraphFracsLabels(dsub['GO_BP'],len(msids))
	if len(x[0])>0:
		printPie(x[0],x[1],'GO:%20Biological%20Processes','GO')

	print '<BR><BR>'
	print "</div>"
	print """<BR><BR><input type=button value = 'Show GO:Molecular Function' name = 'GO_MFDivButton' onclick = 'tog("GO_MFDiv","GO_MFDivButton")'>"""
	print """<div id = 'GO_MFDiv' style="display:none">"""
	f.write(">GO: Molecular Function\n")
	print '<center><h3>GO: Molecular Function</h3></center>'
	subsetGOTable(dsub['GO_MF'],dall['GO_MF'],len(allmsids),len(msids),f,sort="enrichment",linkType="GO",qs=qs,heading="GO_MF",expid=expid,export = export)
	x=getGraphFracsLabels(dsub['GO_MF'],len(msids))
	if len(x[0])>0:
		printPie(x[0],x[1],'GO:%20Molecular%20Function','GO')
	print "</div>"
	f.write(">GO: Cellular Component\n")
	print """<BR><BR><input type=button value = 'Show GO:Cellular Component' name = 'GO_CCDivButton' onclick = 'tog("GO_CCDiv","GO_CCDivButton")'>"""
	print """<div id = 'GO_CCDiv' style="display:none">"""
	print '<center><h3>GO: Cellular Component</h3></center>'
	subsetGOTable(dsub['GO_CC'],dall['GO_CC'],len(allmsids),len(msids),f,sort="enrichment",linkType="GO",qs=qs,heading="GO_CC",expid=expid,export = export)
	x=getGraphFracsLabels(dsub['GO_CC'],len(msids))
	if len(x[0])>0:
		printPie(x[0],x[1],'GO:%20Cellular%20Components','GO')
	print "</div>"
def getMSids(exp_id,c):
	c.execute("""select id from MS where experiment_id=%s""",(exp_id))
	x=c.fetchall()
	return [item[0] for item in x]

		







def printJavaScript():
	print """
	<script type="text/javascript">
	
	function tog(o,button)
	{
	    var e = document.getElementById(o);
	    e.style.display = e.style.display == 'block' ? 'none' : 'block';
	   
	    button = document.getElementsByName(button)[0]
	    if(button.name == 'GOButton')
	    {button.value = button.value == 'Show GO' ? "Hide GO":"Show GO" ;}
	    if(button.name == 'SCANButton')
	    {button.value = button.value == 'Show Scansite' ? "Hide Scansite":"Show Scansite" ;}
	    if(button.name == 'DOMButton')
	    {button.value = button.value == 'Show Domains' ? "Hide Domains":"Show Domains" ;}
	    if(button.name == 'PFAMButton')
	    {button.value = button.value == 'Show PFAM' ? "Hide PFAM":"Show PFAM" ;}
	    if(button.name == 'GO_BPButton'||button.name == 'GO_MFButton' ||button.name == 'GO_CCButton' ||button.name =='PFAMTableButton'||button.name == 'DOMTableButton'||button.name == 'SCANTableButton'||button.name == 'KINTableButton' || button.name == 'PELMTableButton')
	    {button.value = button.value == 'Show Table' ? "Hide Table":"Show Table";}
	    if(button.name == 'GO_BPPieButton'||button.name == 'GO_MFPieButton' ||button.name == 'GO_CCPieButton' || button.name == 'PFAMPieButton'||button.name == 'DOMPieButton' ||button.name == 'SCANPieButton'||button.name == 'KINPieButton' ||button.name == 'PELMPieButton')
	    {button.value = button.value == 'Show Graph' ? "Hide Graph":"Show Graph";}
	    if(button.name == 'GO_BPDivButton')
	    {button.value = button.value == 'Show GO:Biological Process'? "Hide GO:Biological Process":"Show GO:Biological Process"}
	    if(button.name == 'GO_MFDivButton')
	    {button.value = button.value == 'Show GO:Molecular Function'? "Hide GO:Molecular Function":"Show GO:Molecular Function"}
	    if(button.name == 'GO_CCDivButton')
	    {button.value = button.value == 'Show GO:Cellular Component'? "Hide GO:Cellular Component":"Show GO:Cellular Component"}
	   
	    
	    
	}
	</script>
	"""


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
	if (form.tog.value =="show GO terms")
	    {
	    form.tog.value = "hide GO terms";
	    }
        else
	    {
	    form.tog.value = "show GO terms";
	    }
        if (form.tog.value =="show PFAM terms")
	    {
	    form.tog.value = "hide PFAM terms";
	    }
        else
	    {
	    form.tog.value = "show PFAM terms";
	    }
        if (form.tog.value =="show Domain terms")
	    {
	    form.tog.value = "hide Domain terms";
	    }
        else
	    {
	    form.tog.value = "show Scansite terms";
	    }
        if (form.tog.value =="show Scansite terms")
	    {
	    form.tog.value = "hide Scansite terms";
	    }
        else
	    {
	    form.tog.value = "show Scansite terms";
	    }
        
	}
	"""
	print '</script>'

# Go table has id(name)//num in sample//num in exp//enrichment
def sortListByElement(list,index):
	
	tmp = [(item[index],item) for item in list]
	
	tmp.sort()
	
	list = [t[1] for t in tmp]
	return list
		

def subsetGOTable(dSub,dAll,num_exp,num_sample,f,sort="enrichment",heading = "ID",linkType='',qs='',expid='',export = True):
	link = ''
	if linkType == "GO":
		link="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term="
	elif linkType=="domain" or linkType=="pfam":
		link = "http://pfam.sanger.ac.uk/family?entry="
	elif linkType=="kinase":
		pass
	elif linkType == "bind":
		pass
	num = 1
	dataList = []
	print '<center><table><tr bgcolor=blue><th>%s</th><th>Num in all</th></tr>' % heading
	for item in dAll.keys():
		if 'GO' in item:
			name = item+' (' +str(GODict.get(item,'')) +')'
		else: name = item
		dataList.append([name, dAll.get(item,0)])
	data=sortListByElement(dataList,num)
	if len(dataList)>1000: 
		dataList = dataList[-500:]
		print "Only first 500 results shown<br>"
	
	gray = False
	null =[]
	
	for item in data:
		if item[0] == '~~~' or item[0]=='null:':
			item[0]='NULL'
			null = item
		if linkType == "domain" or linkType=="pfam":
			if "Pfam-B" in item[0]:
				link = "http://pfam.sanger.ac.uk//pfamb/"
			else:
				link = "http://pfam.sanger.ac.uk/family?entry="
		if linkType !='bind' and linkType!='kinase' and linkType!=None and link!='':
			item[0]='<a target = "_blank" href="%s">'%(link+item[0].split('(')[0])+item[0]+'</a>'
	if null in data:
		data.remove(null)
		data = data + [null]
	for item in data:
		if '<' in item[0]:
			copyterm = item[0].split('>')[1].split('<')[0]
		else: copyterm = item[0]
		if (item[0].strip()!='NULL' and item!=null) or linkType=='pfam' or linkType == "GO" or linkType == 'bind' or linkType=='kinase' or linkType=='pelm':
			if item == null and linkType=='pfam':
				newqs = qs+'&MFIRST0=%s&MFIRST0_value=%s'%(heading,"~~~")
			else:
				newqs = qs+'&MFIRST0=%s&MFIRST0_value=%s'%(heading,copyterm.split()[0].strip())
				
			url = "%ssubset/search.cgi?%s"%(urlpath,newqs)
			url = url.replace('&','code1')
			url = url.replace('?','code2')
			if gray: print '<tr align = center bgcolor=ccccff><td align= left>'
			else: print '<tr align = center><td align= left>'
			if export:
				if not(gray):
					print '%s</td><td>%s</td></tr>' % (item[0]+"&nbsp;"*3+"<a href='%ssubset/select.cgi?exp_id=%s&amp;url=%s' target='_blank'><img border='0' alt='glass' src='%sglass.png' height='20'></a>"%(urlpath,expid,url,imgPath),item[1])
				else:
					print '%s</td><td>%s</td></tr>' % (item[0]+"&nbsp;"*3+"<a href='%ssubset/select.cgi?exp_id=%s&amp;url=%s' target='_blank'><img border='0' alt='glass' src='%sglass.png' height='20'></a>"%(urlpath,expid,url,imgPath),item[1])
			else:
				print '%s</td><td>%s</td></tr>' % (item[0],item[1])
			gray = not(gray)
			f.write(str(copyterm)+"\t"+str(item[1])+"\n")
	if null!=[] and linkType != 'pfam' and linkType !='GO' and linkType != 'bind' and linkType != 'kinase' and linkType !='pelm':
		if gray:
			print '<tr align = center bgcolor=ccccff><td align= left>'
		else:
			print '<tr align = center><td align= left>'
		print '%s</td><td>%s</td></tr>' % (null[0],null[1])
		f.write(str(null[0])+"\t"+str(null[1])+"\n")
	print '</table></center>'

	


def printPfamDiv(dsub,dall,msids,allmsids,f,qs='',expid='',export = True):
	dsub = dall
	msids = allmsids
	f.write(">PFAM Sites\n")
	print '<BR>'
	subsetGOTable(dsub['pfam_site'],dall['pfam_site'],len(allmsids),len(msids),f,sort="enrichment",linkType="pfam",qs=qs,expid=expid,heading="pfam_site",export = export)
	x=getGraphFracsLabels(dsub['pfam_site'],len(msids))
	if len(x[0])>0:
		printPie(x[0],x[1],'PFAM Domains','pfam')
	print '<BR>'

	
def printDomDiv(dsub,dall,msids,allmsids,f,qs='',expid='',export= True):
	dsub = dall
	msids = allmsids
	f.write(">PFAM Domains\n")
	print "<BR>"
	subsetGOTable(dsub['domain'],dall['domain'],len(allmsids),len(msids),f,sort="enrichment",heading = "domain",linkType="pfam",qs=qs,expid=expid,export= export)
	x=getGraphFracsLabels(dsub['domain'],len(msids))
	if len(x[0])>0:
		printPie(x[0],x[1],'Domains','dom')
	print '<BR>'
	
def printScanDiv(dsub,dall,msids,allmsids,f,qs,expid='',export = True):
	dsub = dall
	msids = allmsids
	if len(dsub['bind'])>0:
		f.write(">Scansite Bind\n")
		print '<center><h3>Scansite Bind</h3></center>'
		subsetGOTable(dsub['bind'],dall['bind'],len(allmsids),len(msids),f,sort="enrichment",heading = "scansite_bind",qs=qs,expid=expid,export = export,linkType="bind")
		x=getGraphFracsLabels(dsub['bind'],len(msids))
		if len(x[0])>0:
			printPie(x[0],x[1],'Scansite Bind','bind')
		print '<BR><BR>'

	if len(dsub['kinase'])>0:
		f.write(">Scansite Kinase\n")
		print '<center><h3>Scansite Kinase</h3></center>'
		subsetGOTable(dsub['kinase'],dall['kinase'],len(allmsids),len(msids),f,sort="enrichment",heading = "scansite_kinase",qs=qs,expid=expid,export = export,linkType='kinase')
		x=getGraphFracsLabels(dsub['kinase'],len(msids))
		if len(x[0])>0:
			printPie(x[0],x[1],'Scansite Kinase','kinase')
		print '<BR><BR>'
	if len(dsub['pelm'])>0:
		badPelm = 1 #put this here to handle empty for now.  Can add again when we want to display pelm 
	#	f.write(">Pelm Kinase\n")
	#	print '<center><h3>Pelm Kinase</h3></center>'
		#subsetGOTable(dsub['pelm'],dall['pelm'],len(allmsids),len(msids),f,sort="enrichment",heading = "pelm_kinase",qs=qs,expid=expid,export = export,linkType='pelm')
		#subsetGOTable(dsub['pelm'],dall['pelm'],len(allmsids),len(msids),f,sort="enrichment",heading = "pelm_kinase",qs=qs,expid=expid,export = export,linkType='foo')
		#x=getGraphFracsLabels(dsub['pelm'],len(msids))
		#if len(x[0])>0:
		#	printPie(x[0],x[1],'Pelm Kinase','pelm')
		#print '<BR><BR>'


def getGraphFracsLabels(dictionary,nmsids):
	if nmsids == 0: return 
	fracs=[]
	labels = []
	for item in dictionary.keys():
		fracs.append(1.0*dictionary[item]/nmsids)
		labels.append(item)
	tuples =  [(fracs[i],labels[i]) for i in range(len(fracs))]
	tuples.sort()
	tuples.reverse()
	return [item[0] for item in tuples],[item[1] for item in tuples]


def printCSSJS():
	print """<style type="text/css">    
    .pg-normal {
        color: black;
        font-weight: normal;
        text-decoration: none;    
        cursor: pointer;    
    }
    .pg-selected {
        color: black;
        font-weight: bold;        
        text-decoration: underline;
        cursor: pointer;
    }
</style>
"""
	print """ <script type="text/javascript">
	
    var pager = new Pager('results', 10);
    pager.init(); 
    pager.showPageNav('pager', 'pageNavPosition'); 
    pager.showPage(1);
</script>"""	

def getParameters(c,msids):
	numpeps=0
	peps = sets.Set()
	for ms in msids:
		query = """select pep_aligned from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id = %s""" %(ms)
		if query[-1]!=' ':
			c.execute(query)
			x=c.fetchall()
		## add check to make sure the phophopep actually fits the criteria
			if x != None:
				x=[item[0] for item in x]
				for item in x:
					query = """select phosphopep.id from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id=%s and pep_aligned='%s'""" %(ms,item)
					c.execute(query)
					x=c.fetchall()
					peps.add(x[0][0])
					numpeps += 1
	ids = len(msids)
	proteins= sets.Set()
	for ms in msids:
		query = """select protein.id from protein join MS on MS.protein_id=protein.id where MS.id=%s"""%(ms)
		if query[-1]!=' ' and query[-1]!= '=':
			c.execute(query)
			x=c.fetchall()
			if x !=None:
				proteins.add(x[0][0])
	
	numProteins = len(proteins)
	return ids, numpeps, numProteins, msids, list(peps),list(proteins)

def sortListByElement(list,index):
	
	tmp = [(item[index],item) for item in list]
	
	tmp.sort()
	
	list = [t[1] for t in tmp]
	list.reverse()
	return list

def getTypeBySpecies(species,expid,c):
	query = """select site_type,count(*) from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id = phosphopep.id join MS on MS.id = MS_phosphopep.MS_id join protein on MS.protein_id=protein.id where experiment_id=%s and species='%s' group by site_type"""%(expid,species.replace("'",""))
	c.execute(query)
	x=c.fetchall()
	x=[(item[0],item[1]) for item in x]
	return dict(x)


def getNumRejected(expid,c):
	if expid ==None: return ''
	query = """select errorLog from experiment where id = %s"""%expid
	c.execute(query)
	x=c.fetchall()
	if len(x)>0:
		logfile = logPath+x[0][0].strip()
		try:
			f = open(logfile,'r')
			text = f.read()
			f.close()
			return text.count("REJECTED")
		except IOError: return 0
	return 0
main()
