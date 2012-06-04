#!/usr/local/bin/python
from pathname import *
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import sys
import cgitb
if displayPythonErrors:
	cgitb.enable()
sys.path.append(path+'includes')
from template import *
from peptideFunctions import *
from tableFunctions import *
NUM_DATA_DIVS=7
from processQueries import *
try:
	f = open(path+'dictionaries/GOdict.dict.pkl','rb')
	GODict = pickle.load(f)
	f.close()
except IOError:
	GODict = {}
sys.path.append(path+'subset/clusters')
from clusters import *
sys.path.append(path+'includes')
sys.path.append(path)
sys.path.append(path+'dictionaries')
from clusterJS import *
from makeDynamicsClusters import *
from searchJS import *
features = {}
from Arith import hypgeomsummore
def main():
	db = MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
	c=db.cursor()
	print '<html><head>'
	print '<title>Search Results</title>'
	printJavaScript()
	pvalJS()
	form = cgi.FieldStorage()
	exp_id = form.getvalue("exp_id",7)
	file_key = form.getvalue("file_key", "")
	user_file_name = form.getvalue("user_file_name", "") # name of the file the user uploaded
	
	stringency = form.getvalue("strin","low")
	if type(stringency)==list and len(stringency)>0: stringency = stringency[0]
	formdata = ''
	for item in form.keys():
		value = form.getvalue(item,'')
		if type(value)==list:
			for thing in form.getvalue(item,[]):
				value=thing
				formdata += '&amp;'+item+'='+value
		if item != 'letter':
			formdata += '&amp;'+item+'='+value

	cluster_file_path = None
	if (file_key != ""):
		cluster_file_path = hackedFileName(exp_id, file_key)
		if (not os.path.exists(cluster_file_path)):
			print "<div class='error'>Could not load cluster file</div>"
			return
	
                ## check if need to do the clustering deal##
		clusterType = form.getvalue("cluster0",'')
		if clusterType == "ALL" and form.getvalue('MFIRST0','') == 'clusters':
			fdata = formdata.replace('cluster0=ALL','')
			set = form.getvalue('MFIRST0_value','').replace('Cluster set: ','')
			cSet = ClusterSet(cluster_file_path)
			names = cSet.clusterSets[set].keys()
			clusterJS(fdata,names)
		        ## handle with javascript
	

	

	#####################################
	##########get search type############
	dataType = form.getvalue("dataType","both")
	
	selected = form.getvalue("selectedids",'')
	
	select = form.getvalue("select","all")
	
	letter=form.getvalue("letter","ALL")
	data_queries= makeDataQueries(form)
	# if no dataqueries, ignore
	
	if data_queries == []:
		dataType = 'metadata'
	
	meta_queries = makeMDataQueries(form)
	#print meta_queries
	# if no metaqueries, ignore
	if meta_queries == [] or meta_queries == [[('MFIRST0', 'None', '')]]:
		dataType= 'data'
	if dataType == "data": meta_queries = []
	######################################

	#######display search################
	printQuery(dataType,data_queries,meta_queries)
	## write query to file so we can get it in the motif results
	queryString = getQuery(dataType,data_queries,meta_queries)
	querynum = random.randint(0,10000000000)
	open("%squery%s"%(outPath,querynum),"w").write(queryString)
	#####################################

	#############get msids###############
	### get msids that meet the query
	
	msids= getData(dataType,data_queries,meta_queries,exp_id,cluster_file_path,c)
	### if we were choosing only from currently selected, filter
	copyselected = selected
	selected = selected.split(':')
	ids = []
	if select !='all' and selected !=['']:
		for item in msids:
			if str(item) in selected: ids.append(item)
		msids=ids
	copymsids = copy.copy(msids) # need a copy so that when we filter by letter, etc. we know what we started with

	### if we have too many, don't allow non-paginated view
	if len(msids) > 100:
		if letter=="ALL":
			letter = "A" 
	lettermsids = []

	### filter by letter
	if letter != 'ALL':
		for item in msids:
			query  = """select acc_gene from MS join protein on MS.protein_id=protein.id where MS.id = %s and acc_gene regexp '^%s'""" %(item,letter)
			c.execute(query)
			x=c.fetchall()
			if len(x) >0: lettermsids.append(item)
	else: lettermsids = msids	
	######################################
	
	#### get experiment parameters##########################
	allmsids = getMSids(exp_id,c)
	
	myMSids,myPeps,myProteins,msids,peps,pros=getParameters(copymsids,meta_queries)
	if myMSids == 0:
		msids = []
	allMSidsNUM,allPepsNUM,allProteinsNUM,allMSids,allPeps,allProteins = getParameters(allmsids,[])
	if select != 'all' and selected !=['']:
		subMSids,subPeps,subProteins,selected,selectedPeps,selectedPro= getParameters(selected,meta_queries)
		
	else: # if using all, background is everything
		subMSids,subPeps,subProteins,selected,selectedPeps,selectedPro=allMSidsNUM,allPepsNUM,allProteinsNUM,allMSids,allPeps,allProteins

	## specificity to be used in calculations
	spec = getBackground(meta_queries,data_queries)
	print "Specificity: ", spec
	#### end variables for specificity level and parameters ############
	
	#################################################
	# print hidden input type with the msids so the parent frame can access them for subset selection
	ms = ''
	if len(msids)>0:
		ms = str(msids[0])
	
	if len(msids)>1:
		for item in msids[1:]:
			ms += ':'+str(item)
	print '<form name=msform action =""><input type=hidden name="msids" value="%s"></form>' % ms
	##################################################

	#################################################
	###print heading with how many things are in the results	
	if select=='all' or selected == ['']:
		if spec == 'MSid':
			print str(myMSids) + ' out of ' + str(subMSids) + ' total MS ids<BR><BR> '
		if spec == 'phosphopep':
			print str(myPeps) +' out of ' +str(subPeps) + ' total phosphopeptides<BR><BR>'
		if spec == 'protein':
			print str(myProteins) +' out of ' + str(subProteins)+' total proteins<br><br>'
	else:
		if spec == 'MSid':
			print str(myMSids) +' ('+str(subMSids)+' in subset considered)' ' out of ' + str(allMSidsNUM) + ' total MS ids in experiment<BR><BR>'
		if spec == 'phosphopep':
			print str(myPeps) +' ('+str(subPeps)+' in subset considered)' ' out of ' + str(allPepsNUM) + ' total phosphopeps in experiment<BR><BR>'
		if spec == 'protein':
			print str(myProteins) +' ('+str(subProteins)+' in subset considered)' ' out of ' + str(allProteinsNUM) + ' total proteins in experiment<BR><BR>'
	################################################
			
	### div to put enriched features###############
	print "<div id='enriched'></div>"
	### along the way, print enriched features to a file handle, enFile, and t the end print that info to this div
	filename = '/tmp/ptmscout/file'+str(random.randint(0,1000000))
	enFile = open(filename,'w')
	
	### end div to put enriched features
	################################################

	################################################
	## make buttons for exporting motif stuff ##
	strallpeps=''
	for item in allPeps:
		strallpeps = strallpeps+str(item) + ':'
	
	strpeps=''
	for item in peps:
		strpeps = strpeps+str(item)+':'
	strallpros = ''
	for item in pros:
		strallpros = strallpros+str(item)+':'
	strpros = ''
	for item in selectedPro:
		strpros = strpros + str(item)+':'
	print """
	<form target="_blank" method="POST" action = "motifs/motifresults.cgi">
	
	<input type=hidden value = '%s' name='background'>
	<input type=hidden value = '%s' name='foreground'>
	<input type=hidden value = '%s' name = 'bgprotein'>
	<input type = hidden value = '%s' name = 'fgprotein'>
	<input type = hidden value = '%s' name = 'querynum'>
	%s&nbsp;&nbsp;<input type = "submit" value = "Calculate Motif Enrichment">
	&nbsp;&nbsp;&nbsp;
	<input type="radio" name = "motif_type" value = "peptide" checked>By peptide
	&nbsp;&nbsp;&nbsp;
	<input type="radio" name = "motif_type" value = "protein">By protein
	&nbsp;&nbsp;&nbsp; <input type = "checkbox" name ="include"> Include sites from all experiments
	<br><br>
	<input type="checkbox" name="emailresults">Email motif results to:
	&nbsp;&nbsp;&nbsp;
	<input type="textarea" name="emailaddress">
	<input type = "hidden" name = "expid" value = "%s">
	</form>
	"""%(strallpeps,strpeps,strpros,strallpros,querynum,helper("Enrichment#Motif_Enrichment"),exp_id)
	print "<BR>"
	selectedms = ''
	for item in allmsids:
		if item != allmsids[-1]:
			selectedms += str(item)+':'
		else: selectedms += str(item)
	#############make fasta files######
	filename1 = "%sexp%s%s.fasta"%(fastaPath,exp_id,random.randint(0,100000000))	
	makeFasta2(c,getMSids(exp_id,c),filename1,meta_queries = [])
	filename2 = "%sexp%s%s.fasta"%(fastaPath,exp_id,random.randint(0,100000000))
	makeFasta2(c,msids,filename2,meta_queries = meta_queries)
	
	########weblogo.cgi will take in number of fasta file
	print "<h3>Aligned sequence profiles</h3>"
	print "<table><tr><td><center>Background</center></td><td><center>Foreground</center></td></tr>"	
	print "<tr><td><img alt='no background peptides' src='%sweblogo.cgi?filename=%s'></td>"%(graphDir,filename1)
	if len(msids)>0:
		print "<td><img alt='no foreground peptides' src='%sweblogo.cgi?filename=%s'></td></tr></table><BR><BR>"%(graphDir,filename2)
	## end motif stuff ##
	#############################################################################
	letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','ALL']
	print '<center>'
	for l in letters[:-1]:
		print """<a href="search.cgi?letter=%s%s">%s</a>""" % (l,formdata,l)
	print """<a href="search.cgi?letter=%s%s">%s</a>""" % ('ALL',formdata,'ALL')
	print '</center><BR>'

	#### check if all msids from same experiment
	### if not, don't do anything with data

	allSameExp = allSameExperiment(msids,c)

	####
	if len(msids)>0:
		if len(lettermsids)>0:
			makeWholeExpTable(exp_id,c,"no search","low",form,["MSid","protein","sequence","site", "pep_aligned"],"",msids=lettermsids,width='85',meta_queries=meta_queries,newtab=True)
			
		print '<BR><BR><BR>'
		### do we have a cached dictionary?
		files = os.listdir(path+"dictionaries/savedDicts")
		exps = []
		for item in files:
			try: exps.append(int(item[4:]))
			except ValueError: pass
		subdicts = {}
		alldicts = {}
		if len(allmsids)<2000:
			subdicts['MSid'] = makeDictFromMSids(msids,c,stringency=stringency)
			alldicts['MSid'] = makeDictFromMSids(allmsids,c,stringency=stringency)
			subdicts['phosphopep'] = makeDictFromPeps(peps,c,stringency=stringency)
			alldicts['phosphopep'] = makeDictFromPeps(selectedPeps,c,stringency=stringency)
			subdicts["protein"] = makeDictFromProteins(pros,c,stringency=stringency)
			alldicts['protein'] = makeDictFromProteins(selectedPro,c,stringency=stringency)
		elif int(exp_id) in exps:
			alldicts['protein'] = pickle.load(open(path+"dictionaries/savedDicts/prodict%s"%exp_id,'r'))
			alldicts['MSid'] = pickle.load(open(path+"dictionaries/savedDicts/msdict%s"%exp_id,'r'))
			alldicts['phosphopep'] = pickle.load(open(path+"dictionaries/savedDicts/pepdict%s"%exp_id,'r'))
			subdicts['MSid'] = makeDictFromMSids(msids,c,stringency=stringency)
			subdicts['phosphopep'] = makeDictFromPeps(peps,c,stringency=stringency)
			subdicts["protein"] = makeDictFromProteins(pros,c,stringency=stringency)
			
			
		else: subdicts['MSid'], subdicts['phosphopep'],subdicts['protein'],alldicts['MSid'],alldicts['phosphopep'],alldicts['protein']={},{},{},{},{},{}
		
		
		print '<br><font size=5><b>GO Data</b></font><BR><BR>'
		printGODocs(c)
		print "<BR>"
		print '<input type=button name = GOButton value = "Show GO" onclick=tog("GO","GOButton")>'
		print '<div id = "GO" style="display:none">'
		printGODiv(subdicts[spec],alldicts[spec],msids,allmsids,subdicts['MSid'],alldicts['MSid'],selected,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=subProteins,totalPeps=subPeps,totalMSids=subMSids,enFile=enFile)
		print '</div>'

		print '<br><font size=5><b>PFAM Data</b></font><BR><BR>'
		print '<input type=button name = PFAMButton value = "Show PFAM" onclick=tog("PFAM","PFAMButton")>'
		print '<div id = "PFAM" style="display:none">'
		printPfamDiv(subdicts['MSid'],alldicts['MSid'],msids,allmsids,enFile=enFile)
		print '</div>'
		print '<br><font size=5><b>Domain Data</b></font><BR><BR>'
		print '<input type=button name = DOMButton value = "Show Domains" onclick=tog("DOM","DOMButton")>'
		print '<div id = "DOM" style="display:none">'
	#	printDomDiv(subdicts['protein'],alldicts['protein'],msids,allmsids,subdicts['MSid'],alldicts['MSid'],enFile=enFile)
		printDomDiv(subdicts[spec],alldicts[spec],msids,allmsids,subdicts['MSid'],alldicts['MSid'],enFile=enFile)
		print '</div>'

		print '<br><font size=5><b>Predictions</b></font><BR><BR>'
		print '<input type=button name = SCANButton value = "Show Predictions" onclick=tog("SCAN","SCANButton")>'
		print '<div id = "SCAN" style="display:none">'
		printScanDiv(subdicts['phosphopep'],alldicts['phosphopep'],len(strallpeps.split(':'))-1,len(strpeps.split(':'))-1,msids,allmsids,subdicts['MSid'],alldicts['MSid'],enFile=enFile)
		print '</div>'

		print '<br><font size=5><b>Data</b></font><BR><BR>'
		print '<input type=button name = DYNButton value = "Show Data" onclick=tog("DYN","DYNButton")>'
		print '<div id = "DYN" style="display:none">'
		
		printDynamicsDiv(msids,exp_id,subdicts['MSid'],alldicts['MSid'],allmsids,selected,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=subProteins,totalPeps=subPeps,totalMSids=subMSids,enFile=enFile)
	
			
		print '</div>'
		print '<br><BR>'


		filename = entropyPath+'exp%s_cluster%s.pkl'%(exp_id,str(random.randint(0,100)))
		try:
			f = open(filename,'wb')
			pickle.dump(features,f)
			f.close()
		except IOError: pass
		print "<input type='hidden' name='featureDict' value='%s'>"%filename.split('/')[-1]
		
		printEnrichedTable(enFile)
		
		enFile.close()
		
		print '</body></html>'
		
	else:
		print "<input type='hidden' name='featureDict' value='tmp'>"
		print 'No Data Found<BR><BR></body></html>'
		
	
	print "<script type='text/javascript'>"
	for i in range(1,30):
		print "if (document.height !=0){"
		print "parent.document.getElementsByName('sub'+%s)[0].style.height=document.height;}"%(str(i))
	print "</script>"
	printFooter()
	
	
	
def printPval(tableName,radioName,textName,buttonName):
	print """%s Mark as significant with p-value less than <input size=10 type=text name="%s" value="0.01">&nbsp;<input type=button name="%s" value="Go" onClick="pValCorrect('%s','%s','%s','%s');pValCorrect('%s','%s','%s','%s')"><br>"""%(helper("Enrichment#Multiple_Hypothesis_Correction"),textName,buttonName,tableName,radioName,textName,"pval"+tableName,tableName,radioName,textName,"pval"+tableName)
	print '<input type =radio name="%s" value= "None" align="left" checked> No multiple hypothesis correction<br>'%(radioName)
	print '<input type =radio name="%s" value= "BF" align="left">Bonferroni correction <b><font color="red"> (in red)</font></b><br>'%(radioName)
	#print '<input type=radio name="%s" value= "RFDR" align="left">Rough False Discovery Rate Correction <b><font color="green"> (in green)</font></b><br>'%(radioName)
	print '<input type=radio name="%s" value= "FDR" align="left">False Discovery Rate Correction <b><font color="orange"> (in orange)</font></b><br>'%(radioName)
	print "<div name='pval%s'></div><br><br>"%tableName
def printGODiv(dsub,dall,msids,allmsids,subMS,allMS,sort="enrichment",selectedids=[],myProteins=0,myPeps=0,myMSids=0,totalProteins=0,totalPeps=0,totalMSids=0,enFile=None):
	if max(dsub['GO_BP'].values()+[1])>1:
		print '<center><h3>GO: Biological Process</h3></center>'
		printPval("GOBPTable","GOBP","BPval","BPButton")
		subsetGOTable(dsub['GO_BP'],dall['GO_BP'],len(allmsids),len(msids),subMS['GO_BP'],allMS['GO_BP'],sort="enrichment",selectedids=selectedids,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=totalProteins,totalPeps=totalPeps,totalMSids=totalMSids,type='protein',name="GOBPTable",linkType="GO",file = enFile, heading="GO BP")
		x=getGraphFracsLabels(dsub['GO_BP'],len(msids))
		if len(x[1])>2:
			printPie(x[0],x[1],'GO%20Biological%20Processes','GO')
		print '<BR><BR>'
	if max(dsub['GO_MF'].values()+[1])>1:
		print '<center><h3>GO: Molecular Function</h3></center>'
		printPval("GOMFTable","GOMF","MFval","MFButton")
		subsetGOTable(dsub['GO_MF'],dall['GO_MF'],len(allmsids),len(msids),subMS['GO_MF'],allMS['GO_MF'],sort="enrichment",selectedids=selectedids,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=totalProteins,totalPeps=totalPeps,totalMSids=totalMSids,type='protein',name="GOMFTable",linkType="GO",file=enFile,heading="GO MF")
		x=getGraphFracsLabels(dsub['GO_MF'],len(msids))
		if len(x[1])>2:
			printPie(x[0],x[1],'GO:%20Molecular%20Function','GO')
	if max(dsub['GO_CC'].values()+[1])>1:
		print '<center><h3>GO: Cellular Component</center></h3>'
		printPval("GOCCTable","GOCC","CCval","CCButton")
		subsetGOTable(dsub['GO_CC'],dall['GO_CC'],len(allmsids),len(msids),subMS['GO_CC'],allMS['GO_CC'],sort="enrichment",selectedids=selectedids,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=totalProteins,totalPeps=totalPeps,totalMSids=totalMSids,type='protein',name="GOCCTable",linkType="GO",file=enFile,heading="GO CC")
		x=getGraphFracsLabels(dsub['GO_CC'],len(msids))
		if len(x[1])>2:
			printPie(x[0],x[1],'GO:%20Cellular%20Components','GO')
	if max(dsub['GO_CC'].values()+[1])<=1 and max(dsub['GO_BP'].values()+[1])<=1 and  max(dsub['GO_MF'].values()+[1])<=1:
		print "Insignificant GO term enrichment"
def getMSids(exp_id,c):
	c.execute("""select id from MS where experiment_id=%s""",(exp_id))
	x=c.fetchall()
	return [item[0] for item in x]


# big table has 3 columns, one for each go aspect
# Go table has id(name)//num in sample//num in exp//enrichment
def sortListByElement(list,index):
	
	tmp = [(item[index],item) for item in list]
	
	tmp.sort()
	
	list = [t[1] for t in tmp]
	return list

def subsetGOTable(dSub,dAll,num_exp,num_sample,subMS,allMS,sort="enrichment",heading = "ID",selectedids=[],pvalue = 1000,name="",type="MSid",myProteins=0,myPeps=0,myMSids=0,totalProteins=0,totalPeps=0,totalMSids=0,linkType=None,file = open('/tmp/ptmscout/test%s.txt'%random.randint(0,10000),'w')):
	subsize=len(selectedids)
	if sort=='enrichment':
		num = 3
	else:
		num = 1
	dataList = []
	print """<center>%s</center><BR>"""%(helper("Enrichment#Enrichment_Table"))
	print '<center><table name="%s"><tr bgcolor=blue><th>%s</th> <th>Num in subset</th> <th>Num in all</th> <th>Enrichment<br>(In experiment)</th></tr>' % (name,heading)
	numtotal,numpicked=0,0
	if type=="MSid":
		numtotal=totalMSids
		numpicked=myMSids
	if type == "phosphopep":
		numtotal=totalPeps
		numpicked=myPeps
	if type=="protein":
		numtotal=totalProteins
		numpicked=myProteins
	for item in dSub.keys():
		numfgMS = subMS.get(item,0)
		numbgMS = allMS.get(item,0)
		
		numinteresting= dAll.get(item,0)
		numtotal = num_exp
		numpicked = num_sample
		numfound = dSub.get(item,0)
		item = str(item)
		if 'GO' in item:
			name = item+' (' +str(GODict.get(item,'')) +')'
		else: name = item
		if numfound > numinteresting:
			numfound = numinteresting
		enrichment = hypgeomsummore(numinteresting,numtotal,numpicked,numfound)
		#enrichment = hgeom(numinteresting,numtotal,numpicked,numfound)
		enrich=enrichment
		
		dataList.append([name, numfound,numinteresting,enrich,0,numfgMS,numbgMS])
	data=sortListByElement(dataList,num)
	## get link
	link = ''
	if linkType == "GO":
		link="http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term="
	elif linkType=="domain" or linkType=="pfam":
		link = "http://pfam.sanger.ac.uk/family?entry="
	elif linkType=="kinase":
		pass
	elif linkType == "bind":
		pass
	gray = False
	pvalue = 0.01
	for item in data:
		if int(item[2])>1 and int(item[1])>1:
			if item[0]=="~~~":
				item0 = "NULL"
			else:
				item0 = item[0]
			if linkType == "domain" or linkType=="pfam":
				if "Pfam-B" in item[0]: link = "http://pfam.sanger.ac.uk//pfamb/"
				else:
					link = "http://pfam.sanger.ac.uk/family?entry="
			if linkType !='bind' and linkType!='kinase' and linkType!=None:
				item[0]='<a target = "_blank" href="%s">'%(link+item[0].split('(')[0])+item0+'</a>'
			enriched = float(item[3])<=pvalue
			
			formattedP = "%.2e" % item[3]
			base, power = formattedP.split('e')
			pval = base+' x ' + '10<sup>%s</sup>' % str(int(power))
			if gray: print '<tr align = center bgcolor=ccccff><td align= left>'
			else: print '<tr align = center><td align= left>'
			if enriched:
				file.write(heading + ': ' +item0+' ('+str(pval)+')'+'\t+\n')
				features[heading +':'+item0]=pval
				print '<b><font color="black">%s</font></b></td><td><b><font color="black">%s</font></b></td><td><b><font color="black">%s</font></b></td><td><b><font color="black">%s</font></b></td></tr>' % (item[0],item[1],item[2],pval)
			else:
				file.write(heading + ': ' +item0+' ('+str(pval)+')'+'\t-\n')
				print '<font color="black">%s</font></td><td><font color="black">%s</font></td><td><font color="black">%s</font></td><td><font color="black">%s</font></td></tr>' % (item[0],item[1],item[2],pval)
			gray = not(gray)
			#if int(item[2])>2:
				#features[heading+': '+item0]=(item[5]*1.0/max(item[6] ,0.001),int(item[6]))
	print '</table></center>'

def printPfamDiv(dsub,dall,msids,allmsids,selectedids=[],enFile=None):
	if max(dsub['pfam_site'].values()+[1])>1:
		print '<center><h3>PFAM Domains  %s</h3></center>'%(helper("Domains#Pfam_Site"))
		printPval("PFAMTable","PFAM","PFAMpval","PFAMButton")
		subsetGOTable(dsub['pfam_site'],dall['pfam_site'],len(allmsids),len(msids),dsub['pfam_site'],dall['pfam_site'],sort="enrichment",heading = "PFAM SITE",selectedids=selectedids,type='phosphopep',name="PFAMTable",linkType="pfam",file=enFile)
		x=getGraphFracsLabels(dsub['pfam_site'],len(msids))
		if len(x[1])>2:
			printPie(x[0],x[1],'PFAM Domains','pfam')
	else:
		print "<BR>Insignificant Pfam site enrichment"

	
def printDomDiv(dsub,dall,msids,allmsids,subMS,allMS,selectedids=[],enFile=None):
	if max(dsub['domain'].values()+[1])>1:
		print '<center><h3>Domains %s</h3></center>'%(helper("Domains"))
		printPval("DOMTable","DOM","DOMpval","DOMButton")
		subsetGOTable(dsub['domain'],dall['domain'],len(allmsids),len(msids),{},{},sort="enrichment",heading = "Domain",selectedids=selectedids,name="DOMTable",linkType='domain',file=enFile)
		x=getGraphFracsLabels(dsub['domain'],len(msids))
		if len(x[1])>2:	
			printPie(x[0],x[1],'Domains','dom')
	else: print "<br>Insignificant domain enrichment"
	
def printScanDiv(dsub,dall,numAllPeps,numMyPeps,msids,allmsids,subMS,allMS,enFile=None):
	if len(dsub.get('bind',[]))>0:
		
		if max(dsub['bind'].values())!=1:
			print '<center><h3>Scansite Bind %s</h3></center>'%(helper("Main_Page#Scansite_Predictions"))
			printPval("BINDTable","BIND","BINDpval","BINDButton")
			subsetGOTable(dsub['bind'],dall['bind'],numAllPeps,numMyPeps,subMS['bind'],allMS['bind'],sort="enrichment",heading = "Scansite Bind",name="BINDTable",linkType="bind",file = enFile)
			x=getGraphFracsLabels(dsub['bind'],len(msids))
			if len(x[1])>2:
				printPie(x[0],x[1],'Scansite Bind','bind')
				

	if len(dsub.get('kinase',[]))>0 and max(dsub.get('kinase',{}).values()+[1])!=1:
		print '<center><h3>Scansite Kinase %s</h3></center>'%(helper("Main_Page#Scansite_Predictions"))
		printPval("KINTable","KIN","KINpval","KINButton")
		subsetGOTable(dsub['kinase'],dall['kinase'],numAllPeps,numMyPeps,subMS['kinase'],allMS['kinase'],sort="enrichment",heading = "Scansite Kinase",name="KINTable",linkType="kinase",file=enFile)
		x=getGraphFracsLabels(dsub['kinase'],len(msids))
		if len(x[1])>2:
			printPie(x[0],x[1],'Scansite Kinase','kinase')
	if max(dsub.get('bind',{}).values()+[1])==1 and max(dsub.get('kinase',{}).values()+[1])==1:
		print "<BR>Insignificant prediction term enrichment"




def getGraphFracsLabels(dictionary,nmsids):
	if nmsids == 0: return 
	fracs=[]
	labels = []
	for item in dictionary.keys():
			fracs.append(1.0*dictionary[item]/nmsids)
			labels.append(item)
	if len(fracs)==0: fracs=[1]
	if len(labels)==0: labels=['null']
	return fracs,labels

def printDynamicsDiv(msids,exp_id,dsub,dall,allmsids,sort="enrichment",selectedids=[],myProteins=0,myPeps=0,myMSids=0,totalProteins=0,totalPeps=0,totalMSids=0,enFile=None):
	copyms = copy.copy(msids)
        chart = ''
	errors = 'true'
	ids=''
	if len(msids)>1:
		ids += str(msids[0])
		for ms in msids[1:]:
			ids += ':'+str(ms)
	elif len(msids)==1:
		ids = msids[0]
	msids = ids
	#### start printing data
	## make javascript function to toggle errorbars
	print '<script type="text/javascript">'
	print 'function redraw() {'
	print 'var d = document.getElementById("graph")'
	print 'var transform = document.graphform.transform.value;'
	print 'var check = document.graphform.errorbars.checked'
	print ' if (check == false)'
	print """{d.innerHTML = "<img  alt='graph unavailable' src='"""+graphDir+"""makeGraph.py?exp_id=%s&msids=%s&errors=%s&transform="+transform+">"} """ % (exp_id,msids,'false')
	
	print 'else'
	print """ {d.innerHTML = "<img  alt='graph unavailable' src='"""+graphDir+"""makeGraph.py?exp_id=%s&msids=%s&errors=%s&transform="+transform+">"}  """ % (exp_id,msids,'true')


	print '}'
	print '</script>'
	
	print '<center>'
	print '<BR><BR><h2>Data</h2>'
	
	# first see if there is actually any data
	c.execute("""select * from data join MS on data.MS_id=MS.id where experiment_id=%s""",(exp_id,))
	x=c.fetchall()
	
	if len(x)>0:
		
		print '<form name = graphform method = GET action = ""><input type = checkbox name = errorbars onChange=redraw() checked>Include Errobars<br>'
		print '<p>Select Data Transform: &nbsp;&nbsp;<select name = "transform" onChange = "redraw()"><option value = None>None</option><option value = "normToMax">Normalize to Maximum</option>'
		print '</select>'
		print '</form>'
		print '<div id = "graph">'
		print '<img alt = "graph unavailable" src="'+graphDir+'makeGraph.py?exp_id=%s&errors=%s&chart=%s&msids=%s">' % (exp_id,errors,chart,msids)
		
		print '<a name = "drawGraph"></a>'
		print '</div>'
		
	else:
		print '<h2>No Data Found</h2>'
	
	
	print '</center>'
	if getDynType(copyms[0],c)=='line':
		print """<center><h3>Dynamic Enrichment %s</h3></center>"""%(helper("Enrichment#Dyanmics"))
		printPval("DYNTable","DYN","DYNpval","DYNButton")
		subsetGOTable(dsub.get('dynamics',{}),dall.get('dynamics',{}),len(allmsids),len(copyms),dsub.get('dynamins',{}),dall.get('dynamics',{}),sort="enrichment",selectedids=selectedids,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=totalProteins,totalPeps=totalPeps,totalMSids=totalMSids,type='MSid',name="DYNTable",file=enFile,heading="dynamics")
	elif getDynType(copyms[0],c)=='bar':
		print """<center><h3>Quantitative Enrichment %s</h3></center>"""%(helper("Enrichment#Dyanmics"))
		printPval("DYNTable","DYN","DYNpval","DYNButton")
		subsetGOTable(dsub.get('dynamics',{}),dall.get('dynamics',{}),len(allmsids),len(copyms),dsub.get('dynamics',{}),dall.get('dynamics',{}),sort="enrichment",selectedids=selectedids,myProteins=myProteins,myPeps=myPeps,myMSids=myMSids,totalProteins=totalProteins,totalPeps=totalPeps,totalMSids=totalMSids,type='MSid',name="DYNTable",file=enFile,heading="dynamics")
		
def specificity(type):
	if type in ['gi','swissprot','entrez_protein','GO_BP','GO_CC','GO_MF','acc_gene','species','domain']:
		return 'protein'
	elif type in ['data','clusters','pfam_site']:
		return 'MSid'
	else:
		return 'phosphopep'

def getBackground(meta_queries,data_queries):
	types=[]
	if data_queries!=[]:
		types.append('MSid')
		return "MSid"
	types=types+[specificity(query[0][1]) for query in meta_queries]
	# use most specific
	if 'phosphopep' in types: return 'phosphopep'
	elif 'MSid' in types: return 'MSid'
	else: return 'protein'


def getParameters(msids,meta_queries):
	#print "<BR><BR>"
	numpeps=0
	peps = []
	msids = list(Set(msids))
	for ms in msids:
		query = """select phosphopep.id from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id = %s""" %(ms)
		if query[-1]!=' ':
			c.execute(query)
			x=c.fetchall()
		## add check to make sure the phophopep actually fits the criteria
			if len(x)>0:
				x=[item[0] for item in x]
				for item in x:
					#print item
					if checkPep2(item,meta_queries,c):
						peps.append(item)
						numpeps += 1
				
		    		
	ids = len(msids)
	proteins= Set()
	for ms in msids:
		query = """select protein.id from protein join MS on MS.protein_id=protein.id where MS.id=%s"""%(ms)
		if query[-1]!=' ' and query[-1]!= '=':
			c.execute(query)
			x=c.fetchall()
			if len(x)>0:
				proteins.add(x[0][0])
	
	numProteins = len(proteins)
	return ids, numpeps, numProteins, msids,list(peps),list(proteins)

def printEnrichedTable(f):
	filename = f.name
	f.close()
	f = open(filename,'r')
	lines = f.readlines()
	lines = [item.split('\t')[0] for item in lines if item.split('\t')[-1].strip()=="+"]
	f.close()
	html= "<b>Features enriched in subset (pvalue < 0.01):</b><br>"
	for item in lines:
		html+= item.strip()
		html+= "<br>"
	if len(lines)==0: html += "None <BR>"
	html += "<br><br>"
	html += "Enriched features are <b>bolded</b> below<BR><BR>"
	print "<script type='text/javascript'>"
	print """document.getElementById('enriched').innerHTML = "%s" """%html
	print "</script>"

def factorial(n):
	#import math
	#return math.sqrt(2*math.pi*n)*(n*1.0/math.e)**n
	if n == 1: return 1
	if n == 0: return 1
	return n*factorial(n-1)

def choose(n,k):
	return factorial(n)/factorial(k)/factorial(n-k)

def hgeom(numinteresting,numtotal,numpicked,numfound):
	try:
		return choose(numinteresting,numfound)*1.0*choose(numtotal-numinteresting,numpicked-numfound)/choose(numtotal,numpicked)*1.0
	except ZeroDivisionError: return 1
	

def makeFasta2(dbConnection,msids,filename,meta_queries = []):
	f = open(filename,'w')
	id = 1
	for item in msids:
		query = """select pep_aligned,phosphopep.id from MS join MS_phosphopep on MS.id = MS_phosphopep.MS_id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id = %s"""%item
		dbConnection.execute(query)
		x=dbConnection.fetchall()
		x = [(thing[0],thing[1]) for thing in x]
		for pep in x:
			seq,id = pep
			if len(fifteenmerize(seq.strip()))==15 and checkPep2(id,meta_queries,dbConnection):
				f.write(">%s\n"%id)
				f.write(fifteenmerize(seq.strip())+"\n")
		id += 1
	f.close()


		
main()
