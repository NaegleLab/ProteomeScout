#!/usr/local/bin/python
from pathname import *
import os
import pylab
import MySQLdb
db = MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
c=db.cursor()
import numpy
import copy
import random
def getMStoGraph2(msids,exp_id,c):
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

			try:
				newname += newnames[-1]
			except IndexError: pass
			
				
		
		x=c.fetchall()
		c.execute("""select phosphopep_id from MS_phosphopep where MS_id=%s""",(ms,))
		try:
			phos = c.fetchall()[0][0]
			c.execute("""select acc_gene from protein join phosphopep on phosphopep.protein_id = protein.id where phosphopep.id = %s""",(phos,))
			y = c.fetchall()
			if len(y)>0: protein = y[0][0]
			else: protein = ''
		except IndexError: protein = ''
		names.append(protein +':'+newname)
	return [(msids[i], names[i]) for i in range(len(msids))]

## input: protein id and experiment ide and database handle
## output: list of tuples of (msid, site). includes multiple site phosphorylation
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
			try:
				newname += newnames[-1]
			except IndexError: pass
			
				
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
	try:
		MS_id=c.fetchall()[0][0]
	except IndexError:
		return [],[]
	c.execute("""select * from data where MS_id=%s""",(MS_id,))
	x=c.fetchall()
	runtypes = [item[2] for item in x]
	# uniquify to get a list of all possible run types
	runtypes = reduce(lambda l, x: x not in l and l.append(x) or l, runtypes, [])
	# get types (same for all runs, so just find for one run type)
	run = runtypes[0]
	# get types
	c.execute("""select type,label from data where MS_id=%s and run = %s order by priority""",(MS_id,run))
	x=c.fetchall()
	types = [(x[i][0],x[i][1]) for i in range(len(x))]
        if len(toAverage(runtypes))>1:
		runtypes.append('Averaged')
	return (runtypes, types)
	
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
        newvec.append(numpy.average([item[i] for item in vectors if item[i]!=None]))
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
			if not l.isalnum() or l.isalpha():
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
			if len(x) == 0: ## if missing data point, CHANGE THIS!
				yvec.append(None)
			else:
				yvec.append(x[0][0])
                runDict['yvector']=yvec

		# if stddevs is empty, then yerror is None, else make a vector of yerrors from stddevs
                yerr = []
                for std in stddevs:
                    c.execute("""select value from data where MS_id=%s and type = %s and label = %s""",(MS_id,std[0],std[1]))
                    x=c.fetchall()
		    try:
			    add = x[0][0]
		    except IndexError: add = 0
                    yerr.append(add)
                yerrs=[]
                for std in newsdevs:
                    if std==0:
                        yerrs.append(0)
                        
                    else:
                        c.execute("""select value from data where MS_id=%s and type = %s and label = %s""",(MS_id,std[0],std[1]))
                        x=c.fetchall()
			try:
				add = x[0][0]
			except IndexError: add = 0
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
            stddevs = stddev(values)
            avDict['yerrors']=stddevs
            dataDict['Averaged']=avDict
            
                
	return dataDict



def plotter(xvec,yvec, yerr = [],color=None, label = ""):
    """ plot line graph, with disconnect for NF values, also add dot for ever non-NF point"""
    ## randomly generate color
    if color == None:
        rgb = (random.random(),random.random(),random.random())
    else: rgb = color
    for i in range(len(yvec)):
	    pylab.errorbar(xvec[i],yvec[i], color = rgb, fmt= "o")
    i = 0
    toPlotX = []
    toPlotY = []
    toPlotYerr = []
    uselabel = label
    ## are there any connected lines?
    allDisconnect = True
    for j in range(len(yvec)-1):
	if yvec[j]!= None and yvec[j+1]!=None:
		allDisconnect = False
    
    while i < len(yvec):
        while yvec[i] != None and i <len(yvec):
            toPlotX.append(xvec[i])
            toPlotY.append(yvec[i])
            if yerr != []:
                toPlotYerr.append(yerr[i])
            i += 1
            if i >=len(yvec): break
        if len(toPlotX)== 1:
		fmt = "o"
		## if there are any connected lines (if any 2 adjacent not None)
		## then don't label the dots
		if not allDisconnect:
			uselabel= ""
		else:
			label = label
        else:
		fmt = "-"
		uselabel= label
	if toPlotX != [] and toPlotY != [] and len(toPlotX)==len(toPlotY):
		if yerr == []:
			pylab.errorbar(toPlotX,toPlotY,yerr = None,color=rgb,fmt = fmt,label=uselabel)
		else:
			pylab.errorbar(toPlotX,toPlotY,yerr = toPlotYerr,color=rgb,fmt = fmt,label=uselabel)
        toPlotX = []
        toPlotY = []
        toPlotYerr = []
        i += 1
    


def dotted(xvec, yvec, yerr = None, color = None, label= ""):
    """ dotted draws dotted line ac if have points abc and b is missing"""
        ## randomly generate color
    if color == None:
        rgb = (random.random(),random.random(),random.random())
    else: rgb = color
    i = 0
    toPlotX = []
    toPlotY = []
    toPlotYerr = []
    dotted = False
    dottedFirst = None
    dottedSecond = None
    while i < len(yvec):
        if not dotted: # in normal line
            if yvec[i] !=None and i < len(yvec)-1:
                toPlotX.append(xvec[i])
                toPlotY.append(yvec[i])
                if yerr != None:
                    toPlotYerr.append(yerr[i])
            else:
                # plot what we have so far
                fmt = "-"
                if len(toPlotX) == 1: fmt = "o"
		if toPlotX != [] and toPlotY != [] and len(toPlotX) == len(toPlotY):
			if yerr == None:
				pylab.errorbar(toPlotX,toPlotY,yerr= None, color = rgb, fmt = fmt, label= label)
			else:
				pylab.errorbar(toPlotX,toPlotY, yerr = toPlotYerr, color = rgb, fmt = fmt, label= label)
                    ## now we are in dotted, set dotted to True and set the first endpoint
                dotted = True
                dottedFirst = max(0,i-1)
            i = i+ 1
        else: # we are in dotted
            ## look for next not None pt.
            fmt = "--"
            while yvec[i] == None:
                i += 1
                if i >= len(yvec):
                    break
            ## found it
            dottedSecond = min(i,len(yvec)-1)
	    if yvec[dottedFirst] == None:
		    yvec[dottedFirst] = 0
	    if yvec[dottedSecond] == None:
		    yvec[dottedSecond] = 0
            ## plot dotted line
            if yerr == None:
                pylab.errorbar([xvec[dottedFirst],xvec[dottedSecond]],[yvec[dottedFirst],yvec[dottedSecond]],yerr= None, color = rgb, fmt = fmt, label= label)
            else:
                pylab.errorbar([xvec[dottedFirst],xvec[dottedSecond]],[yvec[dottedFirst],yvec[dottedSecond]],yerr= [yerr[dottedFirst],yerr[dottedSecond]], color = rgb, fmt = fmt, label= label)
	    if dottedSecond == len(yvec)-1 and yvec[dottedSecond]!=None:
	        fmt = "o"
		if yerr == None:
			pylab.errorbar(xvec[dottedSecond],yvec[dottedSecond],yerr= None, color = rgb, fmt = fmt, label= label)
		else:
			pylab.errorbar(xvec[dottedSecond],yvec[dottedSecond], yerr = yerr[dottedSecond], color = rgb, fmt = fmt, label= label)
            dotted = False
            toPlotX = []
            toPlotY = []
            toPlotYerr = []


def plotBar(xvec, yvec, width, yerr = None, label = "", color = None):
	""" make bar graph, getting rid of NF values"""
	if color == None:
		color = (random.random(),random.random(),random.random())
	newx = []
	newy = []
	newyerr = []
	for i in range(len(yvec)):
		if yvec[i]!=None:
			newx.append(xvec[i])
			newy.append(yvec[i])
			if yerr !=None:
				newyerr.append(yerr[i])
	if yerr == None: newyerr = None
	pylab.bar(newx,newy, yerr = newyerr, label = label, color = color)
	
	
