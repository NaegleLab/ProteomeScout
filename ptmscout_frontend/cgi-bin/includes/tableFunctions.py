from pathname import *
import sys

sys.path.append(path+'includes')
sys.path.append(path+'batch')
from template import *
from peptideFunctions import *
from proteinFcns import *
import copy
import sets

def makeWholeExpTable(exp_id, c,search,stringency,form,cols,url,letter=None,msids=[],width='100',speciesSearch="all",exactProName = None,syns=True,meta_queries=[],icons=True,newtab=False,default = False):
	if url == "sample.cgi": title= "Browse_Experiment"
	elif url == "ambiguity.cgi":title="Ambiguity_Report"
	else: title = "Subset_Selection#Explore_Subset"
	if newtab: tabstring = "target='_blank'"
	else: tabstring = ''
        print "<center>"
    
	# start big table (MS_id, acc, small table with colspan=3)
        if letter != None:
            letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z', 'ALL']
            print '<center>'
            for l in letters:
		print "<a href='%s?exp_id=%s&amp;letter=%s'>%s</a>" % (url,exp_id, l,l)
            print '</center><br>'
        
	print """&nbsp;&nbsp;<a target="_blank" href="%stitle=%s"><img border="0" alt="help" width="20" src="%shelp.jpg"></a>"""%(helpPath,title,imgPath)

	print """<div id="loading">Loading...</div>"""
	
	print '<table  cellpadding="3" cellspacing="0"  width = "95%%" >'
	print '<thead >'
	print '<tr align = "left" bgcolor="#cccccc">'
        if "MSid" in cols: print '<th>MS id</th>'
        if "protein" in cols: print '<th>Protein</th>'
        if "sequence" in cols: print "<th>Trypsinized Peptide</th>"
        
        if "species" in cols: print '<th>Species</th>'
        if "site" in cols: print '<th>Site</th>'
        if "pep_aligned" in cols: print '<th>Aligned Sequence</th>'
        if "pfam" in cols: print '<th>Pfam Site</th>'
        if "predictions" in cols: print '<th>Predictions</th>'
        if "ambiguity" in cols:
		print "<th>Choose protein"
		print """<a href="ambiguity.cgi?expid=%s&default=yes">Select defaults</a>"""%(exp_id)
		if default:
			print """
			<br>
			<font color='red'>*</font> = Default protein assignment different than original.
			"""
		print "</th>"
        print '</tr></thead>'
	if msids!=[]: ms_ids=msids
	elif exactProName !=None:
		query = """select id from protein where acc_gene = '%s'"""%exactProName
		c.execute(query)
		exactPIDs = [item[0] for item in c.fetchall()]
		ms_ids=[]
		for pid in exactPIDs:
			query = """select id from MS where protein_id = %s"""%pid
			c.execute(query)
			x=c.fetchall()
			ms_ids.extend([item[0] for item in x])
		
	elif search == "no search":
		### get total msids, if it is less than 50, display all at once
		query = """select MS.id from MS join protein on MS.protein_id=protein.id where experiment_id=%s"""%exp_id
		c.execute(query)
		x=c.fetchall()
		if len(x)<50:
			pass
	        elif letter != 'ALL':
		    if exp_id!=None:
			    query = """select MS.id from MS join protein on MS.protein_id=protein.id where experiment_id = %s and acc_gene regexp '^%s' order by acc_gene""" % (exp_id,letter)
		    else:
			    query = """select MS.id from MS join protein on MS.protein_id=protein.id where acc_gene regexp '^%s' order by acc_gene""" % (letter)
			    
		else:
		    if exp_id!=None:
			    query = """select MS.id from MS join protein on MS.protein_id=protein.id where experiment_id = %s order by acc_gene"""%(exp_id)
		    else:
			    query = """select MS.id from MS join protein on MS.protein_id=protein.id where order by acc_gene"""
		c.execute(query)
		ms_ids=[item[0] for item in c.fetchall()]
        elif search == "nothing":
		ms_ids=[]
	else:
		if exp_id!=None and exp_id!="None":

			query = "select MS.id from MS join acc on MS.protein_id=acc.protein_id where acc.value regexp '^%s'"%search + " and experiment_id=%s"%exp_id
		else:
			query = "select MS.id from MS join acc on MS.protein_id=acc.protein_id where acc.value regexp '^%s'"%search
			
		c.execute(query)
		x=c.fetchall()
		if len(x)>0:
			tempids = [it[0] for it in x]
		else: tempids = []
		if exp_id!=None and exp_id!="None":
			query = "select MS.id from MS join protein on MS.protein_id=protein.id where acc_gene like '%"+search+"%'" +" and experiment_id=%s"%exp_id
		else:
			query = "select MS.id from MS join protein on MS.protein_id=protein.id where acc_gene like '%"+search+"%'" 
			
		c.execute(query)
		
		x=c.fetchall()
		ms_ids=[]
		if len(x)>0:
			ms_ids = sets.Set([it[0] for it in x] + tempids)
		else: ms_ids=sets.Set(tempids)
		if ms_ids==[]: print "Sorry, your search returned no results."
	
		
	### filter by species
	speciesids = []
	for item in ms_ids:
		query = """select species from protein join MS on MS.protein_id = protein.id where MS.id = %s"""%item
		c.execute(query)
		try:
			species = c.fetchall()[0][0]
		except IndexError: species = ""
		if species == speciesSearch: speciesids.append(item)
	if speciesSearch!="all":
		ms_ids = speciesids

	# for each phosphopep in the experiment, make an item in the table, alternate colors each row
	blue = False
	ambiguous = False
	allms_ids = copy.deepcopy(ms_ids)
	if exp_id == None or exp_id=="None":
		ms_ids = uniqueProteins(ms_ids,c) ### uniquify if for protein search
		if len(ms_ids)>1 or len(ms_ids)==0:
			print "<BR><BR>%s results<BR><BR>"%len(ms_ids)
		else:
			print "<BR><BR>%s result<BR><BR>"%len(ms_ids)
	if len(ms_ids)>500:
		print "Only the first 500 sites have been loaded. Please search by protein to narrow your search"
		ms_ids = ms_ids[0:500]
	ms_ids=list(sets.Set(ms_ids))
	ms_ids = alphabetizeMS(list(ms_ids),c)
	for item in ms_ids:
		query = """select protein_id from MS where id = %s"""%item
		c.execute(query)
		pid=c.fetchall()[0][0]
                sequence, acc, sites,id,species,seq_al,pfam = getTableValues(item, c)
	
		if exp_id==None or exp_id=="None":seq_al,sites = getProteinSites(item,c)
                newseq=[]
                newsites=[]
                if meta_queries!=[]:
			for s in seq_al:
				query = """select phosphopep_id from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id = phosphopep.id where MS_id=%s and pep_aligned = '%s'"""%(item,s)
				c.execute(query)
				pepid = c.fetchall()[0][0]
				if checkPep2(pepid,meta_queries,c):
					newsites.append(sites[seq_al.index(s)])
					newseq.append(s)
			seq_al=newseq
			sites = newsites
        #for site in sites:
                fontTagO=''
                fontTagC=''
                         
                ambiguous = testAmbiguity(sequence,c,species)
                if blue:
                        print '<tr bgcolor="#CCCCCC" valign="top">' 
                else: print '<tr valign = "top" >'
## first column - MS_id
                if "MSid" in cols:
                    print '<td>' #rowspan = %s	>'% str(len(seq_al))
                    print fontTagO
                    print str(item)
                    print fontTagC,'</td>'
                
## second column - print acc_gene column
                if "protein" in cols:
                    print fontTagO	
 
                    if ambiguous:
                        ips = getProteins(sequence,c)
                        proteins = [it[0] for it in ips]
			proteinsIds=orderProteins(pid,item,sequence,c,species=species)
                        ids=[thing[1] for thing in ips]
                        spec=[sp[2] for sp in ips]
                        print """<td onmouseover="showmenu('t%s')" """ % str(item)#acc
                        print """onmouseout="hidemenu('t%s')"> """% str(item)#acc
			#print """rowspan="%s">"""% str(len(seq_al))
                        print """<a %s href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s"><font >"""% (tabstring,urlpath,id,str(exp_id),species)
                        print acc 
                        print '</font></a>'
                        if blue: print '<img src="%sred_flag_blue.jpg" alt="flag">'%imgPath
                        else: print '<img src="%sred_flag.jpg" alt="flat">'%imgPath # will add stuff to this to link out
                        print """<table id='t%s' class="test">""" % str(item) #acc
                        usedaccs=[]
                        for i in range(len(proteins)):
                                # make sure an ambiguity involves the same species
                                if spec[i] ==species and proteins.count(proteins[i])==1 and proteins[i] not in usedaccs:
                                        pro = proteins[i]
                                        proid = ids[i]
					
                                        print """<tr onmouseover='this.style.backgroundColor="\#99ccff";' onmouseout='this.style.backgroundColor="white";'><td><a %s href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">""" % (tabstring,urlpath,proid,str(exp_id),species)
                                        print 	"""%s</a>""" % pro
					if syns:
						for syn in getSynonyms(proid,c):
							print '<br><i>'+syn+'</i>'
					print """</td>"""
                                        usedaccs.append(pro)
                                        # are there multiple isoforms
                                elif spec[i]==species and proteins.count(proteins[i])>1 and proteins[i] not in usedaccs:
                                        want = [prot for prot in proteins if prot==proteins[i]]
                                        count = len(want)
                                        pro = proteins[i]
                                        proid = ids[i]
                                        print """<tr onmouseover='this.style.backgroundColor="\#99ccff";' onmouseout='this.style.backgroundColor="white";'>
                                        <td><a %s href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">""" % (tabstring,urlpath,proid,str(exp_id),species)
                                        print   """%s""" % pro
                                        print """ (%s isoforms)</a>""" % count
					if syns:
						for syn in getSynonyms(proid,c):
							print "<BR><i>"+syn+"</i>"
					print """</td></tr>"""
                                        usedaccs.append(pro)
                                else: pass
                                                                                                                                                        
                        print """<tr bgcolor="\#99ccff" onmouseover='this.style.backgroundColor="white";' onmouseout='this.style.backgroundColor="\99ccff";'><td><a %s href='%scompare.cgi?peptide=%s&amp;exp_id=%s""" % (tabstring,urlpath+'ambiguity/',sequence,str(exp_id))
			print """&amp;species=%s'>View All</a></td></tr>""" % species
                        print "</table>"
                    else:  #not ambiguous
		        if not(syns):
				print '<td >'
				print """<a %s href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">"""% (tabstring,urlpath,id,str(exp_id),species)
				print acc
				print '</a>'
				print "</td>"
			else:
				print """<td onmouseover="showmenu('s%s')" """%item
				
				print """onmouseout="hidemenu('s%s')"> """% str(item)#acc
				print """<a %s href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">"""% (tabstring,urlpath,id,str(exp_id),species)
				print acc
				print '</a>'
				print """<table class="test" id = 's%s'>"""%item
				for syn in getSynonyms(pid,c):
					print '<tr><td>'+'<i>'+syn+'</i></td></tr>'
				print "</table>"
				print '</td>'
## third column starts inset table, each row is a phosphopep (seq_al, site_pos, predictions)
		# first finish row for seq_al[0]
                
                if "sequence" in cols:
                    print '<td>' #rowspan = %s	>'% str(len(seq_al))
                    print fontTagO
                    print sequence
                    print fontTagC,'</td>'
		if "species" in cols:
                    print "<td>%s</td>"%species
                if "site" in cols:
                    print "<td><table>"
                    for i in range(len(sites)):
                        print "<tr>"
                        if blue: print '<td bgcolor="#cccccc">'
                        else: print "<td>"
                        print fontTagO
                        print sites[i]
                        if isKinase(seq_al[i],c) and icons:
				if exp_id!=None and exp_id!="None":
					print "<br>"
					print """<img height=20 src="%skinase.jpg" alt="kinase">"""%imgPath
					if isActive(sites[i],c,item)=='true':
						print """<img height=20 src="%sactive.jpg" alt="active">"""%imgPath
					elif isActive(sites[i],c,item)=='maybe':
						print """<img height=20 src="%squestion.jpg" alt = "maybe active">"""%imgPath
                        print fontTagC
                        print "</td>"
                        print "</tr>"
                    print "</table></td>"
                if "pep_aligned" in cols:
                    print "<td><table>"
                    for i in range(len(seq_al)):
                        print "<tr>"           
                        if blue: print '<td bgcolor="#cccccc">'
			else: print '<td>'
                        print fontTagO
                        print seq_al[i]
                        print fontTagC,'</td>'
                        print "</tr>"
                    print "</table></td>"
                if "pfam" in cols:
                # make pfam column
                    if blue: print '<td bgcolor="#cccccc">'
                    else: print '<td>'
                    print fontTagO
                    if pfam[0] != '~~~':
			print pfam[0]
                    print fontTagC,'</td>'
		if "ambiguity" in cols:
			if blue: print "<td bgcolor='#cccccc'>"
			else: print "<td>"
			print "%s</td>"%makeSelectBox(proteinsIds,c,name='sel%s'%item,default=default)
		if "predictions" in cols:	
		        if blue: print '<td bgcolor="#cccccc">'
                        else: print '<td>'
			
## new table to hold predictions			
			if blue:
				print '<table bgcolor=cccccc width="100%">'
			else: print """<table width="100%">"""
			# get phosphopep id of seq_ali]
			c.execute("""select phosphopep_id from MS_phosphopep join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id=%s and site_pos = %s""",(item,str(sites[i][1:])))
			#print item, sites[i][1:]
			x=c.fetchall()
			if len(x)>0:
				pep = x[0][0]
				c.execute("""select source,value,score from phosphopep_prediction where phosphopep_id = %s order by score""",(pep,))
				x=c.fetchall()
			else: x = []

			# make list of all scan_site kinases containing tuples (value,score)
			scansite_kinases = [(k[1],k[2]) for k in x if k[0]=='scansite_kinase']
			scansite_binds = [(k[1],k[2]) for k in x if k[0]=='scansite_bind' ]
			pelm_kinases = [(k[1],k[2]) for k in x if k[0]=='pelm_kinase' ]
			scansites_other = [(k[1],k[2]) for k in x if k[0]=='scansite']
			# get swissprot to use to link for pelm_kinase
			
			c.execute("""select protein_id from MS where id = %s""",(item,))
			x = c.fetchall()
			if len(x)>0:
				proid = x[0][0]
				c.execute("""select value from acc join phosphopep on acc.protein_id=phosphopep.protein_id where type='swissprot' and acc.protein_id=%s""",(proid,))
				x=c.fetchall()
			else: x = []
			if len(x)>0:
				swissprot=x[0][0]
			else: swissprot=''
			
			def printList(list,name,fontTagO,fontTagC, stringency):
				newlist = copy.deepcopy(list)
				
				if stringency == "High": limit = 0.5
				elif stringency == "Medium": limit = 1.0
				else: limit = 100000000
				touse = []
				for item in newlist:
					if item[0]=='~~~':
						newlist.remove(item)
				if name != 'Pelm Kinase':
					for item in newlist:
						if float(item[1]) <= float(limit):
							touse.append(item)
				else: touse=newlist
				
				if len(touse)>0:
					if name=='Pelm Kinase':
						print '<tr align="left"><th>'+fontTagO+name+fontTagC+'</th></tr>'
						for k in touse:
							if k[1]!=None:
								if swissprot != "":
									print '<tr align="left"><td><a href="http://phospho.elm.eu.org/cgimodel.py?accession='+str(swissprot) +'" target="_blank"'+' >'+fontTagO+k[0]+fontTagC+'</a></td></tr>'
								else:
									print '<tr align="left"><td><a href="http://phospho.elm.eu.org/" target="_blank"'+' >'+fontTagO+k[0]+fontTagC+'</a></td></tr>'
									 
							else:
								print '<tr align="left"><td>'+fontTagO+k[0]+fontTagC+'</td></tr>'
								
								
					else:
						print '<tr align="left"><th>'+fontTagO+name+fontTagC+'</th></tr>'
						for k in touse:
							print '<tr align="left"><td>'+fontTagO+k[0]+' ('+str(k[1])+')'+fontTagC+'</td></tr>'
			printList(scansite_kinases,'Scansite Kinase',fontTagO,fontTagC,stringency)
			printList(scansite_binds, 'Scansite Bind',fontTagO,fontTagC,stringency)
			printList(scansites_other,'Scansite',fontTagO,fontTagC,stringency)
			printList(pelm_kinases,'Pelm Kinase',fontTagO,fontTagC,stringency)
			if len(scansite_kinases)+len(scansite_binds)+len(pelm_kinases)+len(scansites_other)==0:
				
				print '<tr align = "center"><td>-----</td></tr>'
			print '<tr><td></td></tr>'
			print '</table>' # end prediction table
			print '</td>' # end prediction column
			
				

                

		
		
		
		print '</tr>'
		
		blue=not(blue)
	
	
	print '</table><BR><BR>' # end big table
        print "</center>"

def getSynonyms(pid,c):
	## return a list with the protein name and list of synonyms
	query = """select name from protein where id = %s"""%pid
	c.execute(query)
	try:
		x=c.fetchall()[0][0]
		name = x
	except IndexError:
		name = ""
	query = """select value from acc where protein_id = %s and type = '%s'"""%(pid,"gene_synonym")
	c.execute(query)
	x=c.fetchall()
	synonyms = [name]+[item[0] for item in x]
	return synonyms
def makeSelectBox(proteinsIds,c,name="",default = False):
	diffFromDefault = False
	## get default protein id
	## get current protein id
	try:
		defaultid = getDefault(proteinsIds,c)
		current = proteinsIds[0][1][0]
		if defaultid !=current:
			diffFromDefault = True
		if default:
			## put default at top of list instead of current choice
			toRemove = None
			for item in proteinsIds:
				if item[1][0] == defaultid:
					toRemove = item
			if toRemove != None:
				proteinsIds.remove(toRemove)
				proteinsIds = [toRemove] + proteinsIds
		
	except IndexError:
		current = None
	string= "<select name=%s>"%name
	for item in proteinsIds:
		if getAcc(item[1][0],c)!=None:
			string+= "<option value='%s'>%s (%s)</option>"%(getAcc(item[1][0],c),item[0],getAcc(item[1][0],c))
	string+= "</select>"
	if default and diffFromDefault:
		print "<font color='red'>*</font>"
	return string

def orderProteins(pid,ms,sequence,c,species = None):
	pid = getProtein(ms,c)
	ips = getProteins(sequence,c,species=species)
	original=None
	for pro in ips: 
		if pid == pro[1]: original=pro
	if original in ips: ips.remove(original)
	addBack=[]
	for item in ips:
		if item[0]==original[0]: #remove ones with same acc as original, add them back on in beginning
			addBack.append(item)
	for item in addBack: ips.remove(item)
	
	proteins = [it[0] for it in ips]
	proteinsIds = [(proteins[i],ips[i]) for i in range(len(ips))]
	proteinsIds.sort()

	### sort addBack
	accs= []
	newAddBack = []
	for item in addBack:
		accs.append(item[0])
		count = len([True for thing in accs if thing.split('(')[0].strip()==item[0]])
		#print item[0], count, "<br>"
		#count = accs.count(item[0])
		if count == 0: newAddBack.append(item[0])
		else:
			newAddBack.append(item[0]+ ' (%s)' %count)
		
	addBack = [(newAddBack[i],addBack[i]) for i in range(len(addBack))]
	### add back the original
	proteinsIds = [(original[0],original)]+addBack+proteinsIds
	proteins = [im[1][0] for im in proteinsIds]
	newproteins = []
	for prot in proteins:
		count = len([True for thing in newproteins if thing.split('(')[0].strip()==prot])
		if count==0: newproteins.append(prot)
		else: newproteins.append(prot+ ' (%s)'%(count+1))

	#print newproteins

	proteinsIds = [(newproteins[i],proteinsIds[i][1][1:]) for i in range(len(newproteins))]
	return proteinsIds

def getProtein(ms,c):
	query="""select protein_id from MS where id = %s"""%ms
	c.execute(query)
	x=c.fetchall()
	try: 
		return x[0][0]
	except IndexError: 
		return None

def uniqueProteins(msids,c):
	peps = sets.Set()
	pros = sets.Set()
	for item in msids:
		#query = """select pep_aligned from phosphopep join MS_phosphopep on phosphopep.id = MS_phosphopep.phosphopep_id where MS_id = %s"""%item
		query = """select protein_id from MS where id = %s"""%item
		c.execute(query)
		x=c.fetchall()
		if len(x)>0: pros.add(x[0][0])
	newms = []
	for item in pros:
		#query = """select MS_id from phosphopep join MS_phosphopep on phosphopep.id = MS_phosphopep.phosphopep_id where pep_aligned = '%s'"""%item
		query = """select id from MS where protein_id=%s"""%item
		c.execute(query)
		x=c.fetchall()
		if len(x)>0: newms.append(x[0][0])
	
	return newms

def getProteinSites(msid,c):### get all things on that protein... which is what we care about
	query = """select protein_id from MS where id = %s"""%msid
	c.execute(query)
	x=c.fetchall()
	if len(x)>0: protein = x[0][0]
	else:
		protein = None
		return [],[]
	###get pep aligned, site type, site number
	query = """select pep_aligned, site_type, site_pos from phosphopep where protein_id=%s"""%protein
	c.execute(query)
	x=c.fetchall()
	if len(x)>0:
		sites= [item[1]+str(item[2]) for item in x]
		pep_aligned = [item[0] for item in x]
		return pep_aligned,sites
	else: return [],[]
