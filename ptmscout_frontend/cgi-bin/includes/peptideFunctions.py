def makeFasta2(c,msids,filename,meta_queries = []):
	f = open(filename,'w')
	id = 1
	for item in msids:
		query = """select pep_aligned,phosphopep.id from MS join MS_phosphopep on MS.id = MS_phosphopep.MS_id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id = %s"""%item
		c.execute(query)
		x=c.fetchall()
		x = [(thing[0],thing[1]) for thing in x]
		for pep in x:
			seq,id = pep
			if len(fifteenmerize(seq.strip()))==15 and checkPep2(id,meta_queries,c):
				f.write(">%s\n"%id)
				f.write(fifteenmerize(seq.strip())+"\n")
		id += 1
	f.close()
def fifteenmerize(peptide):
	if len(peptide)==15: return peptide
	if 'y' in peptide:
		residue = 'y'
	elif 't' in peptide:
		residue = 't'
	elif 's' in peptide:
		residue = 's'
	elif 'k' in peptide:
		residue = 'k'
	else:
		residue = ''
	if residue != '':
		peptide = peptide.strip()
		middle = peptide.find(residue)
		if middle == 7:
			if len(peptide) <15:
				return peptide+(15-len(peptide))*'.'
			else:
				return peptide[0:15]
		elif middle < 7:
			pep= '.'*(7-middle)+peptide
			if len(pep)>15: return pep[0:15]
			else: return pep + (15-len(pep))*'.'
		else:
			pep = peptide[middle-7:]
			if len(pep)>15: return pep[0:15]
			else:
				return pep + (15-len(pep))*'.'
			
		
	return peptide


def makeFasta(c,expid,filename):
	f = open(filename,'w')
	query = """select pep_aligned from MS join MS_phosphopep on MS.id = MS_phosphopep.MS_id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where experiment_id=%s"""%expid
	c.execute(query)
	x=c.fetchall()
	peps = [item[0] for item in x]
	id = 1
	for pep in peps:
		if len(fifteenmerize(pep.strip())) == 15:
			f.write(">%s\n"%id)
			f.write(fifteenmerize(pep.strip())+"\n")
			id += 1
	f.close()
## need to add case if search is NULL
def checkPep2(phosphopep,queries,c):
	compareType = None
	if True in ['AND' in item[0][0] for item in queries]: compareType="AND"
	if True in ["OR" in item[0][0] for item in queries]: compareType="OR"
	
	for q in queries:
		if True in ["AND" in item[0] for item in q]: compareType="AND"
		if True in ["OR" in item[0][1] for item in q]: compareType="OR"
	for q in queries:
		if q[0][1]=='None' or q[0][2]=='None': return True
		type, val = q[0][1],q[0][2]
		if specificity(type)!= 'phosphopep': return True
		elif type in ['clusters','acc_gene']: return True
		elif type in ['scansite_bind','scansite','scansite_kinase']:
			stringency = q[0][3]
			stringencyDict = {'low':5,'medium':1,'high':0.2}
			if val != "NULL":
				query = """select * from phosphopep_prediction where phosphopep_id=%s and source='%s' and value='%s' and score <= %s"""%(phosphopep,type,val,stringencyDict.get(stringency,1))
				c.execute(query)
				x=c.fetchall()
				if compareType=="AND" or compareType==None:
					if len(x)==0: return False
				elif len(x)>0: return True
			else:
				query = """select * from phosphopep_prediction where phosphopep_id=%s and source='%s' and score <= %s"""%(phosphopep,type,stringencyDict.get(stringency,1))
				c.execute(query)
				x=c.fetchall()
				if compareType=="AND" or compareType==None:
					if len(x)==0: return True
					else: return False
				elif len(x)==0: return True
		elif type in ['pelm_kinase']:
			if val != "NULL":
				query = """select * from phosphopep_prediction where phosphopep_id=%s and source='%s' and value='%s'"""%(phosphopep,type,val)
				c.execute(query)
				x=c.fetchall()
				if compareType=="AND" or compareType==None:
					if len(x)==0: return False
				elif len(x)>0: return True
			else:
				query = """select * from phosphopep_prediction where phosphopep_id=%s and source='%s'"""%(phosphopep,type)
				c.execute(query)
				x=c.fetchall()
				if compareType=="AND" or compareType==None:
					if len(x)==0: return True
					else: return False
				elif len(x)==0: return True
				
				
		elif type in ['pep_aligned','hprd']:
			val = val.replace('Hydrophilic','RNDCEQHKPSTY')
			val = val.replace('Hydrophobic','AFGILMPUVW')
			val = val.replace('Basic','HKR')
			val=val.replace("pY",'y')
			val=val.replace("pT","t")
			val=val.replace("pS","s")
			val=val.replace("X",".")
			query = """select pep_aligned from phosphopep where id=%s and pep_aligned regexp binary '%s'"""%(phosphopep,val)
			c.execute(query)
			x=c.fetchall()
			if compareType == "AND" or compareType==None:
				if len(x)==0: return False
			elif len(x)>0: return True
		else:
			query = """select %s from phosphopep where id = '%s'"""%(type,phosphopep)
			
			c.execute(query)
			x=c.fetchall()
			if x!=None:
				if compareType=="AND" or compareType==None:
					if val not in x[0][0]: return False
				elif compareType=="OR":
					if val in x[0][0]: return True
	if compareType=="OR":
		return False
	return True

def getAcc(pid,c):
	accs=[]
	query="""select value,type from acc where protein_id= %s and type <> 'gene_synonym' and type <> 'genbank'"""%pid
	c.execute(query)
	x=c.fetchall()
	for item in x:
		if item[1].strip() == "swissprot":
			accs.append(item[0])
	for item in x:
		if item[1].strip() == "gi":
			accs.append(item[0])
	for item in x:
		if item[1].strip() == "refseq":
			accs.append(item[0])

	for item in x:
		if item[1].strip() == "ipi":
			accs.append(item[0])
	if len(accs)>0:
		return accs[0]
	else:
		return ""

def specificity(type):
	if type in ['gi','swissprot','entrez_protein','GO_BP','GO_CC','GO_MF','acc_gene','species','domain']:
		return 'protein'
	elif type in ['data','clusters']:
		return 'MSid'
	else:
		return 'phosphopep'

def alphabetizeMS(msids,c):
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


def isKinase(pep_aligned,c):
	query = """select pfam_site from phosphopep where pep_aligned = '%s'""" % (pep_aligned)
	c.execute(query)
	x=c.fetchall()
	x=x[0][0]
	if 'Pkinase' in x: return True
	else: return False
def isActive(site,c,ms_id):
	# get sequence
	import re
	maxSeparation = 35 #maximum separation we will allow for non-exact matches before returning maybe
	query = """select sequence, protein.id from protein join MS on MS.protein_id=protein.id where MS.id = %s""" % ms_id

	c.execute(query)
	result = c.fetchall()
	sequence,protein_id=result[0][0],result[0][1]
	site_num = int(site[1:])
	domainStart,domainEnd = returnKinaseDomainIndexesOfInterest(c,site_num,protein_id)
	# find indices of active region

	DFGindex, APEindex, DFGmatch, APEmatch = findActiveLoopIndexes(sequence, domainStart,domainEnd)
	if DFGindex == -1 or APEindex == -1:
		return 'false'

	if site_num < APEindex and site_num > DFGindex:
		if APEindex - DFGindex < maxSeparation:
			return 'true'
		elif DFGmatch == 'DFG' and APEmatc=='APE':
			return 'true'
		else:
			return 'maybe'
	else:
		return 'false' #if site does not fall within 

def findActiveLoopIndexes(sequence, domainStart, domainEnd): ##need to limit to look at domain start and end
	import re
	# find indices of active region
	DFG = False
	APE = False
	DFGmatch = ''
	APEmatch = ''
	DFGindex = sequence.find('DFG',domainStart,domainEnd)
	APEindex = -1
	minSeparation = 15 #distance, at minimu, we might expect to look downstream of DFG for APE
	domain = sequence[domainStart:domainEnd]
	if DFGindex != -1:
		DFG = True
		DFGmatch = 'DFG'
		APEindex = sequence.find('APE',DFGindex+minSeparation, domainEnd)
		if APEindex != -1:
			APE = True
			APEmatch = 'APE'
	if DFG and APE:
		return DFGindex, APEindex, DFGmatch, APEmatch
	if not DFG:
		p = re.compile('D[FPLYWM]G')
		m = p.findall(domain)
		if m:
			DFGmatch = m[0] #going to assume with high conservation that this will match only one time 
			DFGindex= sequence.find(DFGmatch,domainStart,domainEnd)
		if DFGindex != -1:
			DFG = True
	if not APE and DFG: # only need to look for APE if we found DFG
		p = re.compile('[ASP][PILW][ED]')
		subseq = sequence[DFGindex+minSeparation:domainEnd]
		m = p.findall(subseq)
		if m:
			APEmatch = m[0]
			APEindex = sequence.find(APEmatch,DFGindex+minSeparation,domainEnd) #start looking at DFG then
		else:
			APEindex = -1
	return DFGindex, APEindex, DFGmatch, APEmatch

def returnKinaseDomainIndexesOfInterest(c,pos,protein_id):
	#this is needed for when there are two, or more, possible domains of interest, either because there is a split or physically multiple domains
	start = -1
	end = -1
	query = """select start, stop from protein join domain on domain.protein_id=protein.id where protein.id = %s and domain.label regexp 'kinase' and domain.start < %s and domain.stop > %s""" %(protein_id,pos, pos)

	c.execute(query)
	x= c.fetchall()
	if len(x)>0:
		start, end = x[0][0], x[0][1]
	else:
		start, end = 0,0
	return start, end


# sequence, acc_gene, site_position, site_type,pfam = getTableValues(ms_id,database_cursor)
# Takes an MS id and a database cursor and outputs the phosphopep sequence, acc name of the corresponding protein, position of the sequence, and type (Y, S, or T)
def getTableValues(ms_id, c):
	"""for ms_id, returns a tuple of sequence, acc_gene, site_position, site_type"""
	## get protein_id and sequence
	c.execute("""select protein_id, phosphopep from MS where id = %s""",(ms_id,))
	try:
		x=c.fetchall()
		proid, sequence = x[0][0],x[0][1]
	except IndexError: proid, sequence = None, ''
	## get acc_gene, species
	c.execute("""select acc_gene,species from protein where id=%s""",(proid,))
	try:
		x=c.fetchall()
		acc,species = x[0][0],x[0][1]
	except IndexError: acc, species = '',''
	
	c.execute("""select pep_aligned from phosphopep join MS_phosphopep on phosphopep.id=MS_phosphopep.phosphopep_id where MS_id = %s""",(ms_id,))
	x=c.fetchall()
	seq_al=[item[0] for item in x]
	
	sites = []
	pfam = []
	for seq in seq_al:
		c.execute("""select site_pos,site_type,pfam_site from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where pep_aligned = %s and MS_id=%s""",(seq,str(ms_id)))
		x=c.fetchall()
		sites.append(x[0][1]+str(x[0][0]))
		pfam.append(x[0][2])
	return sequence, acc, sites, proid, species,seq_al,pfam


################ ambiguity functions ##########

# boolean = testAmbiguity(sequence, database_cursor)
# tells whether or not a given sequence might belong to more than one possible protein. Returns true if there is some ambiguity.
def testAmbiguity(sequence,c,species):
	sequence = sequence.upper()
	c.execute("""select protein_id from ambiguity where peptide = %s""",(sequence,))
	x=c.fetchall()
	protein_ids = [item[0] for item in x]
	same_spec=[]
	for id in protein_ids:
		c.execute("""select species from protein where id = %s.""",(id,))
		spec=c.fetchall()[0][0]
		if spec == species:
			same_spec.append(id)
	x = same_spec
	if len(x) > 1:
		return True
	return False

# list of proteins = getProteins(sequence,database_cursor):
# returns a list of possible proteins (given by acc_gene) that might contain a given sequence
def getProteins(sequence,c,species = None):
	sequence = sequence.upper()
	if species == None:
		c.execute("""select acc_gene,protein_id,species from protein join ambiguity on ambiguity.protein_id=protein.id where peptide=%s""",(sequence,))
	else:
		query = """select acc_gene,protein_id,species from protein join ambiguity on ambiguity.protein_id=protein.id where peptide='%s' and species = '%s'"""%(sequence,species)
		c.execute(query)
	x=c.fetchall()
	proteins = [item[0] for item in x]
	ids = [item[1] for item in x]
	sp = [item[2] for item in x]
	## list tuples of (protein, id)
	ips= [(proteins[i], ids[i],sp[i]) for i in range(len(proteins))]
	ips = reduce(lambda l, x: x not in l and l.append(x) or l, ips,[])
	return ips

def getPeps(expid,c):
    """returns a list of all the peptide ids corresponding to experiment expid"""
    query = """select phosphopep.id from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id join MS on MS.id = MS_phosphopep.MS_id where experiment_id = %s""" % expid
    c.execute(query)
    x=c.fetchall()
    peps =[item[0] for item in x]
    return peps

def getAmbiguousPeps(peptides,c,expid):
    """returns a list of all peptide ids from experiment expid that are ambiguous"""
    ambiguousPeps = []
    for peptide in peptides:
        query = """select species,phosphopep from protein join MS on MS.protein_id=protein.id join MS_phosphopep on MS.id = MS_phosphopep.MS_id  where phosphopep_id = %s and experiment_id=%s""" % (peptide,expid)
        c.execute(query)
        x = c.fetchall()
        if len(x)>0:
            species = x[0][0]
            sequence = x[0][1]
            ambiguous = testAmbiguity(sequence, c, species)
            if ambiguous:
                ambiguousPeps.append((peptide,sequence))
    return ambiguousPeps

############################################

def getMSFromPep(peptide, expid,c):
	query="""select MS_id from MS_phosphopep join MS on MS_phosphopep.MS_id=MS.id where experiment_id=%s and phosphopep_id=%s"""%(expid,peptide)
	c.execute(query)
	x=c.fetchall()
	if len(x)>0: return x[0][0]
	else: return None
	
def getMSids(exp_id,c):
	c.execute("""select id from MS where experiment_id=%s""",(exp_id))
	x=c.fetchall()
	return [item[0] for item in x]



def getSiteString(msid,c):
	# get protein name
	query = """select acc_gene from protein join MS on protein.id = MS.protein_id where MS.id = %s"""%msid
	c.execute(query)
	x=c.fetchall()
	try: acc = x[0][0]
	except IndexError: acc = ''
	# get sites
	values = getTableValues(msid,c)
	sites = values[2]
	string = acc
	for item in sites:
		string += "_%s"%item
	return string
