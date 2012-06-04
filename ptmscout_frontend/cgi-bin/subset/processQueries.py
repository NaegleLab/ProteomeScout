from pathname import *
NUM_DATA_DIVS=30
import sys
from sets import Set
sys.path.append(path)
sys.path.append(path+'subset/clusters')

from clusters import *
from template import checkPep2

def printOperator(op):
	print """<span class="operator">%s</span>""" %(op)

def getQuery(dataType,data_queries,meta_queries):
	qstr = ""
	symbols = {'plus':'+','minus':'-','divide':'&divide;','mult':'&times;','geq':'&ge;','leq':'&le;','equals':'=','notequals':'&ne;','lt':'<','gt':'>'}
	stringencyDict = {'low':5,'medium':1,'high':0.2}
	qstr+= "<h3>Foreground search</h3>"
	if dataType not in ["data","both","metadata"]: return 
	if dataType=='data' or dataType == 'both':
		qstr+= "Quantitative Measurement Search<BR>"
		for item in data_queries:
			if "AND" in item[0][0]: printOperator("AND")
			if "OR" in item[0][0]: printOperator("OR")
			for thing in item:
				if thing[1] != 'None':
					if thing[1] not in ['plus','minus','divide','mult','gt','lt','equals','notequals','leq','geq']:
						qstr+= thing[1] + ' '
					else:
						qstr+= symbols.get(thing[1],"")+' '
			qstr+= "<BR>"
		

	if dataType=='both' or dataType=='metadata':
		qstr+= "Metadata Search<BR>"
		meta_queries = [item[0] for item in meta_queries]
		for item in meta_queries:
			if "AND" in item[0]: printOperator("AND")
			if "OR" in item[0]: printOperator("OR")
			if 'scansite' in item[1] and len(item)>3: qstr+= item[1]+': '+str(item[2])+ ' (%s stringency < %s)'%(item[3],stringencyDict.get(item[3],1))+"<BR>"
			else: qstr+= item[1]+': '+str(item[2]) + '<BR>'
		qstr+= "<BR>"
	qstr+= '<BR>'
	return qstr


def printQuery(dataType,data_queries,meta_queries):
	symbols = {'plus':'+','minus':'-','divide':'&divide;','mult':'&times;','geq':'&ge;','leq':'&le;','equals':'=','notequals':'&ne;','lt':'<','gt':'>'}
	stringencyDict = {'low':5,'medium':1,'high':0.2}
	print "<h3>Search Results</h3>"
	if dataType not in ["data","both","metadata"]: return 
	if dataType=='data' or dataType == 'both':
		print "<fieldset><legend>Quantitative Measurement Search</legend>"
		for item in data_queries:
			if "AND" in item[0][0]: printOperator("AND")
			if "OR" in item[0][0]: printOperator("OR")
			for thing in item:
				if thing[1] != 'None':
					if thing[1] not in ['plus','minus','divide','mult','gt','lt','equals','notequals','leq','geq']:
						print thing[1] + ' '
					else:
						print symbols.get(thing[1],"")+' '
			print "<br>"
		print "</fieldset>"

	if dataType=='both' or dataType=='metadata':
		print "<fieldset><legend>Metadata Search</legend>"
		meta_queries = [item[0] for item in meta_queries]
		for item in meta_queries:
			if "AND" in item[0]: printOperator("AND")
			if "OR" in item[0]: printOperator("OR")
			if 'scansite' in item[1] and len(item)>3: print item[1]+': '+str(item[2])+ ' (%s stringency < %s)'%(item[3],stringencyDict.get(item[3],1))+"<br>"
			else: print item[1]+': '+str(item[2]) + '<br>'
		print "</fieldset>"
	print '<BR>'

	
def makeDataQueries(form):
	data_queries = []
	
	for i in range(NUM_DATA_DIVS):
		
		newlist= []
		if i == 0:
			
			for item in form.keys():
				if 'FIRST' in item and 'M' not in item and 'num' not in item:
					if 'num' in form.getvalue(item,""):
						value = 'num:'+form.getvalue(item+'_num',"")
					else: 
						value = form.getvalue(item,"")
					#print item,value
				 	newlist.append((item, value))
		elif form.getvalue("AND"+str(i),"")=='active':
			for item in form.keys():
				if ('AND'+str(i)) in item and 'num' not in item and 'M' not in item:
					if 'num' in form.getvalue(item,""):
						value = 'num:'+form.getvalue(item+'_num',"")
					else: 
						value = form.getvalue(item,"")
				 	newlist.append((item, value))
		elif form.getvalue("OR"+str(i),"")=='active':
			for item in form.keys():
				if ('OR'+str(i)) in item and 'num' not in item and 'M' not in item:
					if 'num' in form.getvalue(item,""):
						value = 'num:'+form.getvalue(item+'_num',"")
					else: 
						value = form.getvalue(item,"")
					newlist.append((item,value))
		countNull = len([item[1] for item in newlist if item[1] =="None"])
		if newlist != [] and countNull < 3:
			if 'FIRST' in newlist[0][0]:
				data_queries.append(newlist)
			else:
				data_queries.append(newlist[1:])
		elif newlist != []:
			if len([item[1] for item in newlist[3:] if item[1] == "None"])==4 and newlist[1][1] in ['gt','lt','geq','leq']:
				newlist = [newlist[0]]+[(newlist[1][0],'plus')]+[(newlist[2][0],'num:0')]+[(newlist[3][0],newlist[1][1])]+[(newlist[4][0],newlist[2][1])]+[(newlist[5][0],'plus')]+[(newlist[6][0],'num:0')]
				if 'FIRST' in newlist[0][0]:
					data_queries.append(newlist)
				else:
					data_queries.append(newlist[1:])
			elif len([item[1] for item in newlist[4:] if item[1]=="None"])==4 and newlist[2][1] in ['gt','lt','geq','leq']:
				newlist = [newlist[1]]+[(newlist[2][0],'plus')]+[(newlist[3][0],'num:0')]+[(newlist[4][0],newlist[2][1])]+[(newlist[5][0],newlist[3][1])]+[(newlist[6][0],'plus')]+[(newlist[7][0],'num:0')]
				data_queries.append(newlist)
				
	#print data_queries
	return data_queries
def makeDataQueries2(form):
	#print form
	data_queries = []
	#print 'here'
	newlist={}
	count=1
	for item in form.keys():
		
		value = '0'
		new=[]
		if 'FIRST' in item and 'M' not in item and 'num' not in item:
			if 'num' in form.getvalue(item,""):
				value = 'num:'+form.getvalue(item+'_num',"")
			else:
				#print 'there'
				value = form.getvalue(item,"")
			x= (item,value)
			new = x
		newlist[count]=new
		count += 1
		

	newlist=newlist.values()
	newlist=[[item for item in newlist if item !=[]]]
	print newlist
	
	countNull = 6-len(newlist)	

	return newlist
def makeMDataQueries(form):
	meta_queries = []
	for i in range(NUM_DATA_DIVS):
		newlist= []
		if i == 0:
			for item in form.keys():
				if 'MFIRST' in item and 'value' not in item:
					val = ""
					if form.getvalue('MFIRST'+str(i)+'_value',"") != "":
						val = form.getvalue('MFIRST'+str(i)+'_value',"")
					if form.getvalue(item,"") == 'acc_gene' or form.getvalue(item,"")=='pep_aligned':
						val = form.getvalue('pro' + str(i),"")
					if form.getvalue(item,"") == 'hprd':
						val = form.getvalue('seq'+str(i),"")
					if form.getvalue(item,"")=='clusters':
						val = form.getvalue('MFIRST'+str(i)+'_value','')+', Cluster Number: '+form.getvalue('cluster'+str(i),"")
					if "scansite" in form.getvalue(item,"").lower():
						stringency = form.getvalue("stringency","low")
						if type(stringency)==list and len(stringency)>0: stringency = stringency[0]
						newlist.append((item,form.getvalue(item,""),val,stringency))
					else:
						newlist.append((item, form.getvalue(item,""),val))
		elif 'active' in form.getvalue("MAND"+str(i),[]):
			for item in form.keys():
				if ('MAND'+str(i)) in item and 'value' not in item:
					val=""
					if form.getvalue('MAND'+str(i)+'_value',"") != "":
						val = form.getvalue('MAND'+str(i)+'_value','')
					if form.getvalue(item,"")=='acc_gene' or form.getvalue(item,"")=='pep_aligned':
						val = form.getvalue('pro' + str(i),"")
					if val !='':
						newlist.append((item,form.getvalue(item,["",""])[1],val))
		elif 'active' in form.getvalue("MOR"+str(i),[]):
			for item in form.keys():
				if ('MOR'+str(i)) in item and 'value' not in item:
					val=""
					if form.getvalue('MOR'+str(i)+'_value',"") != "":
						val = form.getvalue('MOR'+str(i)+'_value',"")
					if form.getvalue(item,"")=='acc_gene' or form.getvalue(item,"")=='pep_aligned':
						val = form.getvalue('pro' + str(i),"")
					if "acc_gene" in form.getvalue(item,[]) or "pep_aligned" in form.getvalue(item,[]):
						val = form.getvalue('pro'+str(i),"")
					newlist.append((item,form.getvalue(item,["",""])[1],val))			
		if newlist != []:
			if newlist[0][1]!="None":
				meta_queries.append(newlist)
	return meta_queries

def getData(dataType,data_queries,meta_queries,exp_id,cluster_file_path,c):
	meta_queries= [item[0] for item in meta_queries]
	dataids=[]
	mdataids =[]
	numq = len(data_queries)
	numqm = len(meta_queries)
	
	if dataType == 'data' or dataType == 'both':
		# find the number of queries
		if numq == 1: 
			
			dataids =  processDataQuery(data_queries[0],exp_id,c) 
			
		elif numq > 1:
			
			# find out whether it is an OR or AND
			
			if 'AND' in data_queries[1][0][0]: 
				logic = 'AND'
			else: 
				logic = 'OR'
			# get the n ms_id sets
			sets = []
			for qnum in range(numq):
				sets.append(processDataQuery(data_queries[qnum],exp_id,c))
			#set1 = processDataQuery(data_queries[0],exp_id,c)
			
			#set2 = processDataQuery(data_queries[1],exp_id,c)
		
			
			if logic == 'OR':
				dataids = sets[0]
				for i in range(1,len(sets)):
					dataids = OR(dataids,sets[i])
				#dataids= OR(set1,set2)
			else:
				dataids = sets[0]
				for i in range(1,len(sets)):
					dataids = AND(dataids,sets[i])
				#dataids= AND(set1,set2)
	
	if dataType == 'metadata' or dataType == 'both':
		if numqm == 1: 
			mdataids= processMDataQuery(meta_queries[0],exp_id,cluster_file_path,c) # for now print, later return it somewhere
			
		elif numqm > 1:
			# find out whether it is an OR or AND
			
			if 'AND' in meta_queries[1][0]: 
				logic = 'AND'
			else: 
				logic = 'OR'
			
			# get the n ms_id sets
			sets = []
			for qnum in range(numqm):
				sets.append(processMDataQuery(meta_queries[qnum],exp_id,cluster_file_path,c))
			#set1 = processMDataQuery(meta_queries[0],exp_id,c)
			
			#set2 = processMDataQuery(meta_queries[1],exp_id,c)
			
			
			if logic == 'OR':
				mdataids = sets[0]
				for i in range(1,len(sets)):
					mdataids = OR(mdataids,sets[i])
				#mdataids= OR(set1,set2)
			else:
				mdataids = sets[0]
				for i in range(1,len(sets)):
					mdataids = AND(mdataids, sets[i])
				#mdataids= AND(set1,set2,metaqueries = meta_queries,c=c)
	if dataType == 'data':
		return dataids
	elif dataType == 'metadata':
		return mdataids
	else:
		return AND(mdataids,dataids)

	
def OR(a,b):
	list = a+b
	list = reduce(lambda l, x: x not in l and l.append(x) or l, list,[])
	return list

def AND(a,b,metaqueries=[],c=''):
	list = []
	for item in a:
		if item in b:
			list.append(item)
	list = reduce(lambda l, x: x not in l and l.append(x) or l, list,[])
	## for each msid, make sure it meets queries
	if metaqueries !=[] and c !='':
		metaqueries=[[it] for it in metaqueries]
		newids = []
		for item in list:
			## get peptide_id
			query="""select phosphopep_id from MS_phosphopep where MS_id=%s"""%item
			c.execute(query)
			x=c.fetchall()
			pepids=[thing[0] for thing in x]
			for pepid in pepids:
				if checkPep2(pepid,metaqueries,c):
					newids.append(item)
		return newids
	return list

def getValue(value,ms,exp_id,c):
	if 'num' in value: 
		return float(value.split(':')[1])
	if value =='None': return 0
	else:
		type, run, label = value.split(':')
		c.execute("""select value from data join MS on data.MS_id=MS.id where MS_id = %s and type = %s and run = %s and label = %s and experiment_id=%s""",(ms,type,run,label,exp_id))
		x=c.fetchall()
		if len(x)>0:
			if x[0][0]==None: return 0
			return float(x[0][0])
	

def operate(value1, op, value2):
	if op == 'plus':
		if value1==None: value1=0
		if value2==None: value2=0
		return value1 + value2
	if op == 'mult':
		if value1==None: value1=0
		if value2==None: value2=0
		return value1 * value2
	if op == 'divide':
		if value1==None: value1=0
		if value2==None: value2=1
		if value2 != 0: return value1*1.0/value2
		else: return 0
	if op == 'minus':
		if value1==None: value1=0
		if value2==None: value2=0
		return value1-value2
	else: return value1

def compare(value1, comp, value2):
	if comp == 'equals': return value1 == value2
	if comp == 'notequals': return value1 != value2
	if comp == 'gt': return value1 > value2
	if comp == 'lt': return value1 < value2
	if comp == 'geq': return value1 >= value2
	if comp == 'leq': return value1 <= value2
	else: return false	
def processDataQuery(query, exp_id,c):
	
	q= [item[1] for item in query]
	for i in range(len(q)):
		if q[i]==None:
			q[i]=0
	
	val1=q[0]
	op1 = q[1]
	val2= q[2]
	comp=q[3]
	val3 =q[4]
	op2 =q[5]
	val4 = q[6]
	query="""select MS.id from MS join protein on MS.protein_id=protein_id where experiment_id = %s order by acc_gene"""%(exp_id)

        
	c.execute("""select MS.id from MS join protein on MS.protein_id=protein.id where experiment_id = %s order by acc_gene""",(exp_id))
	x = c.fetchall()
	
	msids = [item[0] for item in x]
	
	fit = []
	for ms in msids:
		val1m = getValue(val1,ms,exp_id,c)
		val2m = getValue(val2,ms,exp_id,c)
		val3m = getValue(val3,ms,exp_id,c)
		val4m = getValue(val4,ms,exp_id,c)
		if compare(operate(val1m,op1,val2m),comp,operate(val3m,op2,val4m)):
			fit.append(ms)
	return fit

def processMDataQuery(mquery,exp_id,cluster_file_path,c):
	stringencyDict = {'low':5,'medium':1,'high':0.2}
	type=mquery[1]
	value = mquery[2]
	c.execute("""select MS.id from MS join protein on MS.protein_id=protein.id where experiment_id = %s order by acc_gene""" ,(exp_id))
	x = c.fetchall()
	
	msids = [item[0] for item in x]
	query = ' '
	if value !="NULL":
		if type in ['clusters']:
			set = value.split(',')[0][13:]
			number = value.split(',')[1].replace(' Cluster Number: ','')
			if number == "ALL":
				print "Sorry, the view all option for viewing clusters is not yet implemented.<BR><BR><BR>"
			else:
				clusterSet = ClusterSet(cluster_file_path).clusterSets.get(set,{})
				return clusterSet.get(number,[])

		if type in ['gi','swissprot','entrez_protein']:
			query = """select MS.id from MS join acc on MS.protein_id =acc.protein_id where type='%s' and value = '%s' and experiment_id=%s"""%(type,value,exp_id)
		if type in ['GO_BP','GO_MF','GO_CC']:
			query = """select MS.id from MS join protein_GO on MS.protein_id=protein_GO.protein_id join GO on protein_GO.GO_id=GO.id where experiment_id=%s and GO = '%s' """%(exp_id,value)

		if type in ['site_type','pfam_site']:
			query = """select MS.id from MS join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id = phosphopep.id where %s = '%s' and experiment_id=%s"""%(type,value,exp_id)

		if type in ['scansite_bind','scansite_kinase','scansite']:
			stringency = mquery[3]
			query = """select MS.id from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id join phosphopep_prediction on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where source='%s' and experiment_id=%s and value = '%s'  and score <= %s"""%(type,exp_id,value,stringencyDict.get(stringency,1))

		if type in ['pelm_kinase']:
			query = """select MS.id from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id join phosphopep_prediction on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where source='pelm_kinase' and experiment_id=%s and value = '%s'"""%(exp_id,value)
		if type in ['acc_gene']:
			query = "select MS.id from MS join protein on MS.protein_id=protein.id where acc_gene like '%"+value+"%'" +" and experiment_id=%s"%exp_id
		if type in ['species']:
			query = """select MS.id from MS join protein on MS.protein_id = protein.id where experiment_id=%s and %s='%s'"""%(exp_id,type,value)

		if type in ['domain']:
			query = """select MS.id from MS join domain on MS.protein_id=domain.protein_id where experiment_id=%s and label = '%s' """%(exp_id,value)
		if type in ['pep_aligned','hprd']:
			value = value.replace('Hydrophilic','RNDCEQHKPSTY')
			value = value.replace('Hydrophobic','AFGILMPUVW')
			value = value.replace('Basic','HKR')
			value=value.replace("pY",'y')
			value=value.replace("pT","t")
			value=value.replace("pS","s")
			value=value.replace("X",".")
			testval = ''
			for i in range(len(value)):
				if value[i]=='X': testval += '.'
				else: testval += value[i]
			query = """select MS.id from MS join MS_phosphopep on MS.id = MS_phosphopep.MS_id join phosphopep on MS_phosphopep.phosphopep_id = phosphopep.id where experiment_id=%s and pep_aligned regexp binary '%s' """%(exp_id,testval)
		c.execute(query)
		x=c.fetchall()
		if len(x)>0: return [item[0] for item in x]
		else: return []
	else: # value = "NULL"
		query = """select id from MS where experiment_id=%s"""%exp_id
		c.execute(query)
		x=c.fetchall()
		msids = Set([item[0] for item in x])
		if type in ["GO_BP","GO_CC","GO_MF"]:
			query = """select MS.id from GO join protein_GO on GO.id=protein_GO.GO_id join MS on MS.protein_id = protein_GO.protein_id where aspect="P" and experiment_id=%s group by MS.id;"""%exp_id
			c.execute(query)
			x=c.fetchall()
			GOBP = Set([item[0] for item in x])
			query = """select MS.id from GO join protein_GO on GO.id=protein_GO.GO_id join MS on MS.protein_id = protein_GO.protein_id where aspect="F" and experiment_id=%s group by MS.id;"""%exp_id
			c.execute(query)
			x=c.fetchall()
			GOMF = Set([item[0] for item in x])
			query = """select MS.id from GO join protein_GO on GO.id=protein_GO.GO_id join MS on MS.protein_id = protein_GO.protein_id where aspect="C" and experiment_id=%s group by MS.id;"""%exp_id
			c.execute(query)
			x=c.fetchall()
			GOCC = Set([item[0] for item in x])
			
			## get sets of msids with either a GO_BP, GO_CC, or GO_MF, then use set difference 
			if type == "GO_BP":
				msids = msids.difference(GOBP)
			if type == "GO_MF":
				msids = msids.difference(GOMF)
			if type == "GO_CC":
				msids = msids.difference(GOCC)
		if type in ['scansite_bind','scansite_kinase','pelm_kinase']:
			if type == "scansite_bind":
				stringency = mquery[3]
				query = """select MS.id from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id = MS_phosphopep.phosphopep_id join MS on MS.id = MS_phosphopep.MS_id where source = 'scansite_bind' and score <= %s and experiment_id=%s group by MS.id"""%(stringencyDict.get(stringency,1),exp_id)
				c.execute(query)
				x=c.fetchall()
				bind = Set([item[0] for item in x])
				## now check that any peptides from the MSids we DID choose have null pelm
				toRemove = []
				for msid in bind:
					query = """select phosphopep_id from MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id  where MS.id =%s"""%msid
					c.execute(query)
					x=c.fetchall()
					for item in x:
						peptide = item[0]
						query = """select * from phosphopep_prediction where source = 'scansite_bind' and phosphopep_id=%s and score <=%s"""%(peptide, stringencyDict.get(stringency,1))
						c.execute(query)
						x=c.fetchall()
						if len(x) == 0:
							toRemove.append(msid)
				for msid in toRemove:
					if msid in bind:
						bind.remove(msid)
				msids = msids.difference(bind)
			if type == "scansite_kinase":
				stringency = mquery[3]
				query = """select MS.id from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id = MS_phosphopep.phosphopep_id join MS on MS.id = MS_phosphopep.MS_id where source = 'scansite_kinase' and score <=%s and experiment_id=%s group by MS.id"""%(stringencyDict.get(stringency,1),exp_id)
				c.execute(query)
				x=c.fetchall()
				kinase = Set([item[0] for item in x])
				## now check that any peptides from the MSids we DID choose have null kinase
				toRemove = []
				for msid in kinase:
					query = """select phosphopep_id from MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id  where MS.id =%s"""%msid
					c.execute(query)
					x=c.fetchall()
					for item in x:
						peptide = item[0]
						query = """select * from phosphopep_prediction where source = 'scansite_kinase' and phosphopep_id=%s and score <=%s"""%(peptide,stringencyDict.get(stringency,1))
						c.execute(query)
						x=c.fetchall()
						if len(x) == 0:
							toRemove.append(msid)
				for msid in toRemove:
					if msid in kinase:
						kinase.remove(msid)
				msids = msids.difference(kinase)
			if type == "scansite":
				stringency = mquery[3]
				query = """select MS.id from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id = MS_phosphopep.phosphopep_id join MS on MS.id = MS_phosphopep.MS_id where source = 'scansite' and score <=%s and experiment_id=%s group by MS.id"""%(stringencyDict.get(stringency,1),exp_id)
				c.execute(query)
				x=c.fetchall()
				scan = Set([item[0] for item in x])
				## now check that any peptides from the MSids we DID choose have null scansite
				toRemove = []
				for msid in scan:
					query = """select phosphopep_id from MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id  where MS.id =%s"""%msid
					c.execute(query)
					x=c.fetchall()
					for item in x:
						peptide = item[0]
						query = """select * from phosphopep_prediction where source = 'scansite' and phosphopep_id=%s and score <=%s"""%(peptide,stringencyDict.get(stringency,1))
						c.execute(query)
						x=c.fetchall()
						if len(x) == 0:
							toRemove.append(msid)
				for msid in toRemove:
					if msid in scan:
						scan.remove(msid)
				msids = msids.difference(scan)
			if type == "pelm_kinase":
				query = """select MS.id from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id = MS_phosphopep.phosphopep_id join MS on MS.id = MS_phosphopep.MS_id where source = 'pelm_kinase' and experiment_id=%s group by MS.id"""%exp_id
				c.execute(query)
				x=c.fetchall()
				pelm = Set([item[0] for item in x])
				## now check that any peptides from the MSids we DID choose have null pelm
				toRemove = []
				for msid in pelm:
					query = """select phosphopep_id from MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id  where MS.id =%s"""%msid
					c.execute(query)
					x=c.fetchall()
					for item in x:
						peptide = item[0]
						query = """select * from phosphopep_prediction where source = 'pelm_kinase' and phosphopep_id=%s"""%item
						c.execute(query)
						x=c.fetchall()
						if len(x) == 0:
							toRemove.append(msid)
				for msid in toRemove:
					if msid in pelm:
						pelm.remove(msid)
				msids = msids.difference(pelm)
		return list(msids)
			
