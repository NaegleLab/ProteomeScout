# pid: protein id
# species
# type: null of mouse, for human 'tissue' or 'cell'
# sizex,sizey size of graph

# make global column lists
from pathname import *
import MySQLdb
db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
c=db.cursor()

c.execute("""describe expression_human""")
x=c.fetchall()
EXPRESSION_HUMAN_COLS=[item[0] for item in x][:-1]

c.execute("""describe expression_mouse""")
x=c.fetchall()
EXPRESSION_MOUSE_COLS = [item[0] for item in x][:-1]

c.execute("""describe expression_NCI60""")
x=c.fetchall()
EXPRESSION_NCI60_COLS = [item[0] for item in x][:-1]

def getProbeListFromPeptide(peptide,c):
	pass
def getProbeOptions(protein_id,c):
	c.execute("""select probeset_id from protein_expression where protein_id=%s""",(protein_id,))
	x=c.fetchall()
	try:
		probes=[item[0] for item in x]
	except IndexError: probes = []
	return probes

def getDefaultProbe(protein_id,c):
	# might make this more sophisticated later
	getProbeOptions(protein_id,c)[0]

## when graphing multipy, probeid is a list of probeids
def makeGraph(probeid, species, c,table='tissue',columns = 'all',multiple=False):
	print '<SCRIPT type="text/javascript" LANGUAGE="JavaScript1.2">'
	if not(multiple):
		if table=="cell":
			print 'document.write("<h2>%s<\/h2><h4>By Cell Type<\/h4>")' % probeid
		else:
			print 'document.write("<h2>%s<\/h2><h4>Tissue Expression Levels<\/h4>")' % probeid
	else:
		if table=="cell":
			print 'document.write("<h4>By Cell Type<\/h4>")' 
		else:
			print 'document.write("<h4>Tissue Expression Levels<\/h4>")' 
										
		
	print 'g = new BAR_GRAPH("hBar");'
	# mouse
	if species == 'mus musculus':
		if columns == 'all':
			indices = [ i for i in range(1,len(EXPRESSION_MOUSE_COLS))]
		else:
			indices=[]
			for col in columns:
				indices.append(EXPRESSION_MOUSE_COLS.index(col))
		if not(multiple):
			c.execute("""select * from expression_mouse where probeset_id=%s""",(probeid,))
			data=c.fetchall()[0]
		 ### begin multiple probes stuff
		else:
			datasets=[] # datasets is a list of lists of values, one list for each probe
			for probe in probeid:
				c.execute("""select * from expression_mouse where probeset_id=%s""",(probe,))
				x=c.fetchall()
				if len(x)!=0:
					data = x[0]
					vals=[]
					for item in indices:
						vals.append(data[item])
					datasets.append(vals)
				else:
					vals = []
					for item in indices:
						vals.append(0)
					datasets.append(vals)
		### end multiple probes stuff                   
				
		values= []
		for item in indices:
			values.append(data[item])
		labels = ''
		for item in indices[:-1]:
			labels = labels  + EXPRESSION_MOUSE_COLS[item] + ','
		if len(indices)>0: labels = labels  +EXPRESSION_MOUSE_COLS[indices[-1]]
	elif table=='tissue':
		# get data from table
		# c.execute("""select * from expression_human where probeset_id=%s""",(probeid,))
		# data=c.fetchall()[0]
		if columns == 'all' or columns == ["all"]:
			indices = [ i for i in range(1,len(EXPRESSION_HUMAN_COLS))]
		else:
			indices=[]
			#print 'document.write("%s");'% columns
			for col in columns:
				try:
					indices.append(EXPRESSION_HUMAN_COLS.index(col))
				except ValueError: pass
		if not(multiple):
			c.execute("""select * from expression_human where probeset_id=%s""",(probeid,))
			data=c.fetchall()[0]
			values= []
			for item in indices:
				values.append(data[item])
		### begin multiple probes stuff
		else:
			datasets=[] # datasets is a list of lists of values, one list for each probe
			for probe in probeid:
				c.execute("""select * from expression_human where probeset_id=%s""",(probe,))
				x=c.fetchall()
				if len(x)!=0:
					data = x[0]
					vals=[]
					for item in indices:
						vals.append(data[item])
					datasets.append(vals)
				else:
					vals = []
					for item in indices:
						vals.append(0)
					datasets.append(vals)
		### end multiple probes stuff
				
		labels = ''
		for item in indices[:-1]:
			labels = labels  + EXPRESSION_HUMAN_COLS[item] + ','
		if len(indices)>0: labels = labels  +EXPRESSION_HUMAN_COLS[indices[-1]]		
	elif table == 'cell':
# get data from table
		#c.execute("""select * from expression_NCI60 where probeset_id=%s""",(probeid,))
		#2Bdata=c.fetchall()[0]
		if columns == 'all' or columns == ["all"]:
			indices = [ i for i in range(1,len(EXPRESSION_NCI60_COLS))]
		else:
			indices=[]
			for col in columns:
				try:
					indices.append(EXPRESSION_NCI60_COLS.index(col))
				except ValueError: pass

		if not(multiple):
			c.execute("""select * from expression_NCI60 where probeset_id=%s""",(probeid,))
			values = []
			try:
				data=c.fetchall()[0]
				for item in indices:
					values.append(data[item])
			except IndexError:
				print 'document.write("No data found");'
				
		 ### begin multiple probes stuff
		else:
			datasets=[] # datasets is a list of lists of values, one list for each probe
			for probe in probeid:
				c.execute("""select * from expression_NCI60 where probeset_id=%s""",(probe,))
				x=c.fetchall()
				if len(x)!=0:
					data = x[0]
					vals=[]
					for item in indices:
						vals.append(data[item])
					datasets.append(vals)
				else:
					vals = []
					for item in indices:
						vals.append(0)
					datasets.append(vals)
					### end multiple probes stuff         
		labels = ''
		for item in indices[:-1]:
			labels = labels  + EXPRESSION_NCI60_COLS[item] + ','
		if len(indices)>0: labels = labels  +EXPRESSION_NCI60_COLS[indices[-1]]
	else:
		values = []
		labels = ''
	# make values into a string:
	if not(multiple):
		valuesStr =''
		for v in values[:-1]:
			valuesStr = valuesStr + str(v) + ','
		if len(values)>0: valuesStr = valuesStr + str(values[-1])
	elif datasets != []:
		print 'g.barWidth=3;'
		valuesStr=''
		# is is the index of the column number
		for i in range(len(datasets[0])):
			# j is the index of the probe number
			val=''
			for j in range(len(probeid)):
				if j==0:
					#a=str(datasets[j][i])
					#b=a.split('.')[0]#+'.'+a.split('.')[1][0]
					val += str(datasets[j][i])
					if j == len(probeid)-1:
						val+= ','
						
				elif j==len(probeid)-1:
					#a=str(datasets[j][i])
					#b=a.split('.')[0]#+'.'+a.split('.')[1][0]
					val += ';'+str(datasets[j][i]) + ','
				else:
					#2Ba=str(datasets[j][i])
					#b=a.split('.')[0]#+'.'+a.split('.')[1][0]
					val+= ';'+str(datasets[j][i])
			valuesStr+=val
		valuesStr=valuesStr[:-1]
	else:
		valuesStr=''
		
			
	
	print 'g.values = "'+valuesStr+'";'
	print 'g.labels="' + labels + '";'
	print 'g.barWidth=15;'
	if columns=='all' or len(indices)>30:
		print 'g.charts=2';
	elif len(indices)>20:
		print 'g.charts=2';
	print 'g.labelBorder="0px";'
	print 'g.labelBGColor="";'
	print 'g.percValuesSize="0";'
	if multiple:
		print 'g.showValues=1;'
		print 'g.barWidth=10;'
		print 'g.labelSpace=5;'
		legendStr=''
		for item in probeid:
			c.execute("""select acc_gene from protein join protein_expression on protein.id=protein_expression.protein_id where probeset_id=%s""",(item,))
			x=c.fetchall()[0][0]
			legendStr += item + ' ('+x+')'
			if probeid.index(item) != len(probeid)-1: legendStr+=','
		print 'g.legend = "'+legendStr+'";'
		print 'g.titleSpace = 10;'
	#print 'g.barLength=2;'
	print 'g.showValues=1'
	print 'g.titles="Cell Type, ,Expression";'
	print 'g.titleColor="white";'
	print 'g.titleBGColor="blue";'
	print 'g.labelAlign="left";'
	print 'g.titleBorder = "0px";'
	#if valuesStr!='':
	print 'document.write(g.create());'
	print '</SCRIPT>'
