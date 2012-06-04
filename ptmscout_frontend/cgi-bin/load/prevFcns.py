from sets import Set
def allZeros(list):
	if len([item for item in list if item == 9]) == len(list): return True
	else: return False
def getType(text):
	"""getType(text) returns ([run labels],[types]) where types are drom the data:type:label column headings"""
	dcols,rcol = getColumns(text)
	types,labels = getTypes(dcols,rcol,text)
	if rcol!=None:
		runs = Set([item.split('\t')[rcol] for item in text.split('\n')if item !=''])
	else: runs = ['average']
	return (list(runs),[item[0] for item in types])
def getTestLines(peptide,text):
	"""getTestLines(peptide,text), gets all the lines with a given peptide, the lines we will use to test the graphing function"""
	lines = text.split('\n')
	testLines = []
	for item in lines:
		if peptide in item: testLines.append(item)
	return testLines
def pepColumn(text):
	"""pepColumn(text) returns the index of the column containing peptides (pep:tryps or peps: heading)"""
	cols = text.split('\n')[0]
	cols= [item.strip() for item in cols.split('\t')]
	for item in cols:
		if "pep" in item: return cols.index(item)
def previewPeptide(text,pepCol):
	"""previewPeptide(text,pepCol) takes text and the index of the peptide column and grabs the first peptide it finds to test"""
	line = text.split('\n')[2]
	return line.split('\t')[pepCol]
def getColumns(text):
	"""getColumns(text) returns a list of columns containing data and the index of the run column in the form ([datacols],runcol)"""
	if text == '':
		return None,None
	cols = text.split('\n')[0]
	cols= [item.strip() for item in cols.split('\t')]
	dataColumns = []
	runCol = None
	for item in cols: 
		if 'data' in item.lower(): dataColumns.append((item,cols.index(item)))
		if 'run' in item.lower(): runCol = cols.index(item)
	return dataColumns, runCol

def getTypes(dataColumns,runCol,text):
	"""getTypes(dataColumns,runCol,text) returns ([types],[labels])"""
	types = [(item[0].split(':')[1],item[1]) for item in dataColumns]
	# uniquify to get a list of all possible run types
	types = reduce(lambda l, x: x not in l and l.append(x) or l, types, [])
	labels = [(item[0].split(':')[2],item[1]) for item in dataColumns]
	# uniquify to get a list of all possible run types
	labels = reduce(lambda l, x: x not in l and l.append(x) or l, labels, [])
	if runCol==None: runs = [('average',None)]
	else:
		runs = [(item.split('\t')[runCol].strip(),runCol) for item in text.split('\n') if item != '']
		runs = reduce(lambda l, x: x not in l and l.append(x) or l, runs, [])
	if len(toAverage([item[0] for item in runs]))>1 and 'average' not in [item[0] for item in runs]:
		runs.append(('Averaged',None))
	return (types,labels)
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
        newvec.append(avg([item[i] for item in vectors]))
    return newvec
def avg(list):
	if len(list)==0: return 0
	return sum(list)*1.0/len(list)

def std(list):
	differences = [(item-avg(list))**2 for item in list]
	differences = sum(differences)*1.0/len(list)
	return differences**0.5
	
def stddev(vectors):
    for item in vectors:
	    if allZeros(item): vectors.remove(item)
    newvec = []
    for i in range(len(vectors[0])):
        vec= [item[i] for item in vectors]
        newvec.append(std(vec))
    return newvec
        
# get data to graph dictionary of key(run):value(data list of the form below)
# data for a single run: need dictionary with the following keys(xvector, yvector, yerrors, xlabels, title)
# so this is dictionary of dictionaries: key(run): value (dictionary of data stuff) with the keys above
def getDataToGraph(pepLines,types,labels,runCol,dcols):
	dataDict = {}
	try:
		runs = [item.split('\t')[runCol] for item in pepLines]
	except TypeError: runs = []
	for line in pepLines:
		runDict = {}
		try:
			runDict['title']=line.split('\t')[runCol]
		except TypeError:
			runDict['title']='average'
		datatypes = [item for item in types if 'stddev' not in item[0]]   
		stddevs =[item for item in types if 'stddev' in item[0]]
		stddevs = [(item[0].split(':')[2],item[1]) for item in dcols if 'dev' in item[0]]
		labels = [item[0].split(':')[2] for item in dcols if 'dev' not in item[0]]
		# xlabels are the label types
		runDict['xlabels']=labels
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

		
		yvec = []
		for data in types:
			if 'stddev' not in data[0]:
				try:
					yvec.append(float(line.split('\t')[data[1]]))
				except ValueError: yvec.append(0)
			
                runDict['yvector']=yvec
		yerr = []
		for item in stddevs:
			yerr.append(line.split('\t')[item[1]])
                if len(yerr)!= len(yvec):
                    yerr=[0 for item in yvec]
                runDict['yerrors']=yerr

		dataDict[runDict['title']]=runDict
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

		
def getColumns(text):
	if text == '': return None
	cols = text.split('\n')[0]
	cols= [item.strip() for item in cols.split('\t')]
	dataColumns = []
	runCol = None
	for item in cols: 
		if 'data' in item.lower(): dataColumns.append((item,cols.index(item)))
		if 'run' in item.lower(): runCol = cols.index(item)
	return dataColumns, runCol
