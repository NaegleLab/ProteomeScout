from makeDynamicsClusters import *
import sets

def getDomainPep(peps,c):
    # get the proteins
    proteins = []
    for pep in peps:
        c.execute("""select protein_id from phosphopep where id = %s""",(pep))
        x=c.fetchall()
	if len(x)>0:
		proteins.append(x[0][0])
    domains = []
    domDict={}
    for item in proteins:
        c.execute("""select label from domain where protein_id = %s""",(item))
        x=c.fetchall()
	if len(x)>0:
		x=[item[0] for item in x]
		# remove duplicates for a single protein
		x = reduce(lambda l, y: y not in l and l.append(y) or l, x,[])
        domains = domains + list(x)
    for s in domains:
        domDict[s] = domains.count(s)
    return domDict
    
                
def getDomainMS(msids,c):
    """returns a dictionary of {domain: number of proteins that have that domain}"""
    domains = []
    domDict={}
    for item in msids:
        c.execute("""select label from MS join domain on MS.protein_id=domain.protein_id where MS.id = %s""",(item))
        x=c.fetchall()
	if len(x)>0:
		x=[item[0] for item in x]
		# remove duplicates for a single protein
		x = reduce(lambda l, y: y not in l and l.append(y) or l, x,[])
        domains = domains + list(x)
    for s in domains:
        domDict[s] = domains.count(s)
    return domDict
def getDomainProtein(proteins,c):
    """returns a dictionary of {domain: number of proteins that have that domain}"""
    domains = []
    domDict={}
    for item in proteins:
        c.execute("""select label from domain where protein_id = %s""",(item))
        x=c.fetchall()
	if len(x)>0:
		x=[item[0] for item in x]
		# remove duplicates for a single protein
		x = reduce(lambda l, y: y not in l and l.append(y) or l, x,[])
        domains = domains + list(x)
    for s in domains:
        domDict[s] = domains.count(s)
    return domDict

def getGOPep(peps,c):
    # get the proteins
    proteins = []
    for pep in peps:
        c.execute("""select protein_id from phosphopep where id = %s""",(pep))
        x=c.fetchall()
	if len(x)>0:
		proteins.append(x[0][0])
    
    BP, MF, CC = [],[],[]
    for item in proteins:
        countB,countM,countC = 0,0,0
        query = """select GO, aspect,term from GO join protein_GO on GO.id=protein_GO.GO_id where protein_GO.protein_id=%s""" % item
        c.execute(query)
        # make list of tuples (GO id, aspect)
        x=c.fetchall()
        GOs = [(term[0],term[1],term[2]) for term in x]
        for go in GOs:
            if go[1]=='P':
		    BP.append(go[0])
		    countB += 1
            elif go[1]=='F':
		    MF.append(go[0])
		    countM += 1
            else:
		    CC.append(go[0])
	            countC += 1
	if countB == 0: BP.append('null:')
	if countM == 0: MF.append('null:')
	if countC == 0: CC.append('null:')
    GO_BPDict, GO_MFDict, GO_CCDict = {},{},{}
    for item in BP:
        GO_BPDict[item]=BP.count(item)
    for item in MF:
        GO_MFDict[item]=MF.count(item)
    for item in CC:
        GO_CCDict[item]=CC.count(item)
    return GO_BPDict, GO_MFDict, GO_CCDict
def getGOProtein(proteins,c):
	# get the proteins
    BP, MF, CC = [],[],[]
    for item in proteins:
        countB,countM,countC = 0,0,0
        query = """select GO, aspect,term from GO join protein_GO on GO.id=protein_GO.GO_id where protein_GO.protein_id=%s""" % item
        c.execute(query)
        # make list of tuples (GO id, aspect)
        x=c.fetchall()
        GOs = [(term[0],term[1],term[2]) for term in x]
        for go in GOs:
            if go[1]=='P':
		    BP.append(go[0])
		    countB += 1
            elif go[1]=='F':
		    MF.append(go[0])
		    countM += 1
            else:
		    CC.append(go[0])
	            countC += 1
	if countB == 0: BP.append('null:')
	if countM == 0: MF.append('null:')
	if countC == 0: CC.append('null:')
    GO_BPDict, GO_MFDict, GO_CCDict = {},{},{}
    for item in BP:
        GO_BPDict[item]=BP.count(item)
    for item in MF:
        GO_MFDict[item]=MF.count(item)
    for item in CC:
        GO_CCDict[item]=CC.count(item)
    return GO_BPDict, GO_MFDict, GO_CCDict
def getGOD(msids,c):
	# get the proteins
    proteins = []
    for ms in msids:
        c.execute("""select protein_id from MS where id = %s""",(ms))
        x=c.fetchall()
	if len(x)>0:
		proteins.append(x[0][0])
      
    BP, MF, CC = [],[],[]
    for item in proteins:
        countB,countM,countC = 0,0,0
        query = """select GO, aspect,term from GO join protein_GO on GO.id=protein_GO.GO_id where protein_GO.protein_id=%s""" % item
        c.execute(query)
        # make list of tuples (GO id, aspect)
        x=c.fetchall()
        GOs = [(term[0],term[1],term[2]) for term in x]
        for go in GOs:
            if go[1]=='P':
		    BP.append(go[0])
		    countB += 1
            elif go[1]=='F':
		    MF.append(go[0])
		    countM += 1
            else:
		    CC.append(go[0])
	            countC += 1
	if countB == 0: BP.append('null:')
	if countM == 0: MF.append('null:')
	if countC == 0: CC.append('null:')
    GO_BPDict, GO_MFDict, GO_CCDict = {},{},{}
    for item in BP:
        GO_BPDict[item]=BP.count(item)
    for item in MF:
        GO_MFDict[item]=MF.count(item)
    for item in CC:
        GO_CCDict[item]=CC.count(item)
    return GO_BPDict, GO_MFDict, GO_CCDict

def getScanPep(peps,c,stringency="medium"):
    stringencyDict = {'low':5,'medium':1,'high':0.2}
    scan,bind, kinase,pelm = [],[],[],[]
    for pep in peps:
            countscan,countbind,countkin,countpelm = 0,0,0,0
	    c.execute("""select source,value from phosphopep_prediction where phosphopep_id = %s and score <= %s""",(pep,stringencyDict.get(stringency,1)))
	    x=c.fetchall()
	    scans = []
	    if len(x)>0:
		    scans = [(item[0],item[1]) for item in x]
	    for s in scans:
                    item = s
		    if item[0]=='scansite':
			    scan.append(item[1])
                            countscan += 1
		    elif item[0]=='scansite_bind':
			    bind.append(item[1])
                            countbind += 1
		    elif item[0] == 'scansite_kinase':
			    kinase.append(item[1])
                            countkin += 1
		    elif item[0] == 'pelm':
			    pelm.append(item[1])
                            countpelm += 1
            if countscan == 0: scan.append("null:")
            if countbind == 0: bind.append("null:")
            if countkin == 0:
                kinase.append("null:")
            query="""select  value from phosphopep_prediction where phosphopep_id=%s and source='pelm_kinase'"""%pep
            c.execute(query)
            x=c.fetchall()
            if len(x)>0: pelm =pelm+[item[0] for item in x]
            else:
                pelm.append("null:")
    scanDict, bindDict, kinaseDict,pelmDict = {},{},{},{}
    for item in scan: scanDict[item]=scan.count(item)
    for item in bind: bindDict[item]=bind.count(item)
    for item in kinase: kinaseDict[item]=kinase.count(item)
    for item in pelm: pelmDict[item] = pelm.count(item)
    return scanDict, bindDict, kinaseDict,pelmDict

def getScanMS(msids,c,stringency="medium"):
    stringencyDict = {'low':5,'medium':1,'high':0.2}
    scan,bind, kinase,pelm = [],[],[],[]
    for ms in msids:
	    c.execute("""select source,value from phosphopep_prediction join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where MS_phosphopep.MS_id=%s and score < %s""",(ms,stringencyDict.get(stringency,1)))
	    x=c.fetchall()
	    scans = []
	    if len(x)>0:
		    scans = [(item[0],item[1]) for item in x]
	    for s in scans:
                    item = s
		    if item[0]=='scansite':
			    scan.append(item[1])
		    elif item[0]=='scansite_bind':
			    bind.append(item[1])
		    elif item[0] == 'scansite_kinase':
			    kinase.append(item[1])
		    elif item[0] == 'pelm':
			    pelm.append(item[1])
            query="""select  value from phosphopep_prediction  join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where MS_phosphopep.MS_id=%s and source='pelm_kinase'"""%ms
            c.execute(query)
            x=c.fetchall()
            if len(x)>0: pelm =pelm+[item[0] for item in x]
    scanDict, bindDict, kinaseDict,pelmDict = {},{},{},{}
    for item in scan: scanDict[item]=scan.count(item)
    for item in bind: bindDict[item]=bind.count(item)
    for item in kinase: kinaseDict[item]=kinase.count(item)
    for item in pelm: pelmDict[item] = pelm.count(item)
    return scanDict, bindDict, kinaseDict,pelmDict



def makeDictFromProteins(proteins,c,stringency="medium"):
    d= dict()
    d['GO_BP'],d['GO_MF'],d['GO_CC']=getGOProtein(proteins,c)
    d['domain']=getDomainProtein(proteins,c)
    return d

def makeDictFromPeps(peps,c,stringency="medium"):
    d= dict()
    d['scan'],d['bind'],d['kinase'],d['pelm']=getScanPep(peps,c,stringency=stringency)
    d['GO_BP'],d['GO_MF'],d['GO_CC']=getGOPep(peps,c)
    d['domain'] = getDomainPep(peps,c)
    return d    

def makeDictFromMSids(msids,c,stringency="medium",data=True):
    d = dict()
    d['GO_BP'],d['GO_MF'],d['GO_CC']=getGOD(msids,c)
    d['pfam_site']= getPfamSite_MS(0,c,msids=msids)
    if data:
        d['dynamics']=getDataD(msids,c)
    d['domain']=getDomainMS(msids,c)
    d['scan'],d['bind'],d['kinase'],d['pelm'] = getScanMS(msids,c,stringency=stringency)
    return d

def getDataD(msids,c):
    mostUp = []
    mostDown = []
    peak = []
    fold=[]
    early = []
    data = {}
    keys = [] 
    for ms in msids: 
        keys.extend(makeDynamicsFeatures(ms,c).keys())
    data = {}
    for key in keys: data[key]=[]
    for ms in msids:
        dict=  makeDynamicsFeatures(ms,c)
        for item in dict.keys():
            data[item].append(dict[item])
    newdata = {}
    for key in keys:
        for item in data[key]:
            newdata[key+':'+str(item)] = data[key].count(item)
    return newdata

    
def getPfamSite_MS(expid,c,msids=[]):
    if msids == []:
        msids = getMSids(expid,c)
    pfams = []
    for item in msids:
        query ="""select pfam_site from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id = phosphopep.id join MS on MS.id=MS_phosphopep.MS_id where MS.id = %s""" %(item)
        c.execute(query)
        x=c.fetchall()
        if len(x)>0: pfams.append(x[0][0])
    d = {}
    for item in pfams:
        d[item]=pfams.count(item)
    return d

def getPfamDomains_PRO(expid,c):
    domains = []
    ### get domains
    query = """select label from domain join MS on MS.protein_id=domain.protein_id where experiment_id=%s"""%expid
    c.execute(query)
    x=c.fetchall()
    domains = sets.Set([item[0] for item in x])
    d={}
    for item in domains:
        query = """select count(*) from domain join MS on MS.protein_id=domain.protein_id where label='%s' and experiment_id=%s group by domain.protein_id"""%(item,expid)
        c.execute(query)
        x=c.fetchall()
        for thing in x: d[item]=len(x)

    return d
def makeSummaryDict(expid,c):
    d = {}
    d['pfam_site'] = getPfamSite_MS(expid,c)
    d['domain'] = getPfamDomains_PRO(expid,c)
    d['scan'],d['bind'],d['kinase'],d['pelm']=getScanPep(getPeptides(expid,c),c)
    d['GO_BP'],d['GO_MF'],d['GO_CC']=getGOProtein(getProteins(expid,c),c)
    return d
def getMSids(exp_id,c):
	c.execute("""select id from MS where experiment_id=%s""",(exp_id))
	x=c.fetchall()
	return [item[0] for item in x]
def getProteins(expid,c):
    query = """select protein_id from MS where experiment_id=%s"""%expid
    c.execute(query)
    x=c.fetchall()
    return list(sets.Set([item[0] for item in x]))
def getPeptides(expid,c):
    query = """select phosphopep_id from MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id where experiment_id = %s"""%expid
    c.execute(query)
    x=c.fetchall()
    return [item[0] for item in x]
