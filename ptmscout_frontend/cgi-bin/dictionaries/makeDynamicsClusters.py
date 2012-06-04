from pathname import *
import MySQLdb

db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
c=db.cursor()

import sys
import sets

import copy

def getDynamicsVec(MSid,c):
    query = """select experiment_id from MS where id = %s"""%MSid
    c.execute(query)
    expid=c.fetchall()[0][0]
    data = getDataToGraph(MSid,getType(expid),c)
    return data['average']['yvector']

def makeDynamicsFeatures(MSid,c):
    query = """select experiment_id from MS where id = %s"""%MSid
    c.execute(query)
    expid=c.fetchall()[0][0]
    dataList = {}
    runs = getDataToGraph(MSid,getType(expid),c).keys()
    for item in runs:
        dataList[item]=getDataToGraph(MSid,getType(expid),c).get(item,{})
    d = {}
    for run in runs:
        data = dataList.get(run,{})
        vector = data.get('yvector',[])
        xvec =data.get('xvector',[])
        xlabels = data.get('xlabels',[])
        if xvec == None: xvec = range(len(vector))
        if len(vector)>0 and len([item for item in vector if item==None])==0:
            if getDynType(MSid,c)=="line": # also check that it is time here
                try:
                    reg = [(vector[i+1]-vector[i])/(xvec[i+1]-xvec[i]) for i in range(len(xvec)-1)]
                except TypeError: reg = [0]
                d['peakUpReg('+str(run)+')']=str(xvec[reg.index(max(reg))])+'-'+str(xvec[reg.index(max(reg))+1])
                d['peakDownReg('+str(run)+')']=str(xvec[reg.index(min(reg))])+'-'+str(xvec[reg.index(min(reg))+1])
                d['peakTime('+str(run)+')']=xvec[(vector).index(max(vector))]
                d['early('+str(run)+')']= d['peakUpReg('+str(run)+')'][0]==0 and d['peakTime - '+str(run)]==xvec[1]
            else:
                ## need all pairs
                reg = []
                indices = []
                for i in range(len(xvec)):
                    for j in range(len(xvec)):
                        if i !=j:
                            reg.append(abs(vector[j]-vector[i]))
                            indices.append([i,j])
                d['peakDifference('+str(run)+')']=str(xlabels[indices[reg.index(max(reg))][0]])+'-'+str(xlabels[indices[reg.index(max(reg))][1]])
                d['lowestDifference('+str(run)+')']=str(xlabels[indices[reg.index(min(reg))][0]])+'-'+str(xlabels[indices[reg.index(min(reg))][1]])
                d['peakValue('+str(run)+')']=xlabels[(vector).index(max(vector))]

            #d['xvec - '+str(run)]=xvec
            #d['yvec - '+str(run)]=vector
            d['fold('+str(run)+')']=int(max(vector))
            
        else: d['peakValue('+str(run)+')'],d['fold('+str(run)+')'],d['peakUpReg('+str(run)+')'],d['peakDownReg('+str(run)+')'],d['early(' + str(run)+')'] = 0,0,0,0,0
    return d


def getDynType(MSid,c):
    query = """select type from data where MS_id = %s"""%MSid
    c.execute(query)
    x=c.fetchall()
    if len(x)>0:
        for item in x:
            if 'time' in item[0]: return 'line'
    return 'bar'


def getMStoGraph2(msids,exp_id,c):
	names =[]
	for ms in msids:
		c.execute("""select site_type,site_pos from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id = %s""",(ms,))
		x=c.fetchall()
		c.execute("""select phosphopep_id from MS_phosphopep where MS_id=%s""",(ms,))
		phos = c.fetchall()[0][0]
		c.execute("""select acc_gene from protein join phosphopep on phosphopep.protein_id = protein.id where phosphopep.id = %s""",(phos,))
		y = c.fetchall()
		if len(y)>0: protein = y[0][0]
		else: protein = ''
		names.append(protein +':'+x[0][0]+str(x[0][1]))
	return [(msids[i], names[i]) for i in range(len(msids))]

def getMStoGraph(pid,exp_id,c):
	c.execute("""select id from MS where protein_id=%s and experiment_id=%s""",(pid,exp_id))
	x=c.fetchall()
	msids = [item[0] for item in x]
	names =[]
	for ms in msids:
		c.execute("""select site_type,site_pos from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id = %s""",(ms,))
		x=c.fetchall()
		newnames = []
		for item in x:
			newnames.append( item[0]+str(item[1]))
		if len(newnames)==1:
			if newnames[0] in names: newname = newnames[0]+'('+str(names.count(newnames[0])+1)+')'
			else: newname = newnames[0]
		else:
			newname = ''
			for name in newnames[:-1]:
				newname += name+'/'
			
			newname += newnames[-1]
			
				
		names.append(newname)
	return [(msids[i], names[i]) for i in range(len(msids))]
# need to know how many runs, return tuple of ([runs],[types]) where types is [(type, label), ...]
def allZeros(list):
    allz=True
    for item in list:
        if item !=0:
            allz=False
    return allz

def getType(exp_id):
	# look at stuff for one ms_id, since all will (theoretically!) be the same
	c.execute("""select id from MS where experiment_id=%s""",(exp_id,))
	MS_id=c.fetchall()[0][0]
	c.execute("""select * from data where MS_id=%s""",(MS_id,))
	x=c.fetchall()
	runtypes = [item[2] for item in x]
	# uniquify to get a list of all possible run types
	runtypes = reduce(lambda l, x: x not in l and l.append(x) or l, runtypes, [])
        if len(runtypes)>0:
            # get types (same for all runs, so just find for one run type)
            run = runtypes[0]
            c.execute("""select type,label from data where MS_id=%s and run = %s order by priority""",(MS_id,run))
            x=c.fetchall()
            if len(x)>0:
                types = [(x[i][0],x[i][1]) for i in range(len(x))]
            else: types = []
            if len(toAverage(runtypes))>1:
		runtypes.append('Averaged')
        else: return [],[]
        return (runtypes, types)
	
           
	# get types
	
# return True if need to average
def toAverage(runs):
        toAverage=[]
        for item in runs:
            if not(item.isalpha()) and item.isalnum():
                toAverage.append(item)
        return toAverage

# vectors is a list of vectors
def average(vectors):
    for item in vectors:
	    if allZeros(item): vectors.remove(item)
    newvec = []
    for i in range(len(vectors[1])):
        newvec.append(sum([item[i] for item in vectors])/len(vectors))
    return newvec


def stddev(vectors):
    for item in vectors:
	    if allZeros(item): vectors.remove(item)
    newvec = []
    for i in range(len(vectors[0])):
        vec= [item[i] for item in vectors]
        newvec.append(numpy.std(vec))
    return newvec
        
# get data to graph dictionary of key(run):value(data list of the form below)
# data for a single run: need dictionary with the following keys(xvector, yvector, yerrors, xlabels, title)
# so this is dictionary of dictionaries: key(run): value (dictionary of data stuff) with the keys above
def getDataToGraph(MS_id,types,c):
	runs = copy.deepcopy(types[0])
        if 'Averaged' in runs: runs.remove('Averaged')
	types = types[1]
	dataDict = {}
	for run in runs:
		runDict = {}
		
		# title is just run type
		runDict['title']=run

		datatypes = [item for item in types if 'stddev' not in item[0]]
                    
		stddevs =[item for item in types if 'stddev' in item[0]]
                newsdevs=[]
                for d in datatypes:
                    havelabel = False
                    for s in stddevs:
                        if s[1]==d[1]:
                            newsdevs.append(s)
                            havelabel=True
                    if not(havelabel):
                        newsdevs.append(0)

		# xlabels are the label types
		labels = [item[1] for item in types]
		labels = reduce(lambda l, x: x not in l and l.append(x) or l, labels,[])
		runDict['xlabels']=labels

		# if labels are numbers, they should go in the x vector, else the x vector is None
		allnums = True
		for l in labels:
                    if not l.isdigit(): #or l.isalpha():
                        allnums = False
                        break
		if allnums:
			xvec = [int(s) for s in labels]
			runDict['xvector']=xvec
		else:
			runDict['xvector']=None

		# yvector should contain all data points
                yvec = []
		for data in datatypes:
			c.execute("""select value from data where MS_id=%s and type = %s and label = %s and run = %s""",(MS_id,data[0],data[1],run))
			x=c.fetchall()
			if len(x) == 0:
				yvec.append(0)
			else:
				yvec.append(x[0][0])
                runDict['yvector']=yvec

		# if stddevs is empty, then yerror is None, else make a vector of yerrors from stddevs
                yerr = []
                for std in stddevs:
                    c.execute("""select value from data where MS_id=%s and type = %s and label = %s""",(MS_id,std[0],std[1]))
                    x=c.fetchall()
                    add = x[0][0]                        
                    yerr.append(add)
                yerrs=[]
                for std in newsdevs:
                    if std==0:
                        yerrs.append(0)
                        
                    else:
                        c.execute("""select value from data where MS_id=%s and type = %s and label = %s""",(MS_id,std[0],std[1]))
                        x=c.fetchall()
                        add = x[0][0]                        
                        yerrs.append(add)
                runDict['yerrors']=yerrs

		dataDict[run]=runDict
        if len(toAverage(runs))>1:
            avDict = {}
            xvector,labels = dataDict[runs[0]]['xvector'],dataDict[runs[0]]['xlabels']
            avDict['title']='Averaged'
            avDict['xvector']=xvector
            avDict['xlabels']=labels
            # get data values
            values = [dataDict[run]['yvector'] for run in toAverage(runs)]
	    
            #calculate averages
            yvector = average(values)
	    
            avDict['yvector']=yvector
            #stddevs = stddev(values)
            #avDict['yerrors']=stddevs
            dataDict['Averaged']=avDict
            
            
                
	return dataDict



