from pathname import *
import pickle
import copy

f = open(path+'dictionaries/colors.dict.pkl','rb')
colors = pickle.load(f)
f.close()

def makeDict(pid,c):
	c.execute("""select start,stop,label from domain where protein_id=%s""",(pid,))
	x=c.fetchall()
	d = {}
	for item in x:
		start, stop,name = item
		if name not in d.keys():
			d[name]={"start":start, "stop":stop,"name":name}
		else:
			key = name + '%s'% len([thing for thing in d.keys() if name in thing])
			d[key]={"start":start, "stop":stop,"name":name}
	return d
def makepDict(pid,c):
	c.execute("""select id,site_pos,site_type from phosphopep where protein_id=%s""",(pid,))
	x=c.fetchall()
	d = {}
	for item in x:
		id,pos,stype = item
		name = stype+str(pos)
		c.execute("""select source,value from phosphopep_prediction where id=%s""",(id,))
		y=c.fetchall()
		ptypes = {}
		for thing in y:
			ptypes["source"]=thing[0]
			ptypes["value"]=thing[1]
		d[id]={"site_pos":pos, "name":name,"ptypes":ptypes,"id":id}
	return d
def makePExpList(pid,c,exp_id=11):
	c.execute("""select phosphopep_id from MS_phosphopep join MS on MS_phosphopep.MS_id=MS.id where experiment_id=%s""",(exp_id,))
	x=c.fetchall()
	x=[item[0] for item in x]
	ids=[]
	for id in x:
		c.execute("""select protein_id from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id where phosphopep_id=%s""",(id,))
		z=c.fetchall()
		y=[item[0] for item in z]
		if pid in y:
			ids.append(id)
	return ids
def makePOthersList(pid,c,exp_id):
	c.execute("""select id from phosphopep where protein_id = %s""",(pid,))
	x=c.fetchall()
	ids = [item[0] for item in x]
	new=[]
	ex=[]
	for item in ids:
		c.execute("""select experiment_id from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id where phosphopep_id=%s""",(item,))
		x=c.fetchall()
		if len(x)>0:
			x = [id[0] for id in x]
		else: x=[]
			#return [],[]
		if str(exp_id) in [str(thing) for thing in x]:
			ex.append(item)
		else: new.append(item)
	return new,ex

# make a phosphosite dictionary containing only the phosphosites found from the specific experiment
def makepExpDict(pid,c,d,exp_id=11):
	c.execute("""select phosphopep_id from MS where experiment_id=%s""",(exp_id,))
	x=c.fetchall()
	x=[item[0] for item in x]
	ids=[]
	d={}
	for id in x:
		c.execute("""select protein_id from MS where phosphopep_id=%s""",(id,))
		y=c.fetchall()[0]
		if y==pid: ids.append(y)
	for item in ids:
	        id,pos,stype = item
	        name = stype+str(pos)
	        c.execute("""select source,value from phosphopep_prediction where id=%s""",(id,))
	        z=c.fetchall()
		ptypes = {}
		for thing in z:
			ptypes["source"]=thing[0]
			ptypes["value"]=thing[1]
		d[name]={"site_pos":pos, "name":name,"ptypes":ptypes}
	return d

def makepOtherDict(exp,all):
	new = copy.deepcopy(all)
	for item in new.keys():
		if item in exp.keys():
			del new[item]
	return new
																		
	
def getLength(pid,c):
	c.execute("""select sequence from protein where id = %s""",(pid,))
	seq=c.fetchall()[0][0]
	return len(seq)
d=dict()
d['test']={"start":50,"stop":150,"name":"test"}
d['again']={"start":210,"stop":250,"name":"again"}
p=dict()
p['test']={"site_pos":75,"name":'Y75'}
class Domain:
	def __init__(self,dictionary,PLENGTH,BAR, HEIGHT):
		self.d=dictionary
		self.PLENGTH = PLENGTH
		self.BAR=BAR
		self.HEIGHT=HEIGHT
	

	def convert(self,xvalue):
		return int(1.0*xvalue/self.PLENGTH*self.BAR)
	def getCoords(self):
		x1,x2 = self.d.get("start"),self.d.get("stop")
		x1,x2 = self.convert(x1),self.convert(x2)
		width = x2-x1
		return (x1,20,width,self.HEIGHT)
	def getTextRect(self):
		text = self.d.get("name")
		coords=self.getCoords()
		return (text,coords[0],0,coords[2],"center")

def order(keys,c):
	positions=[]
	d={}
	new=[]
	for key in keys:
		c.execute("""select site_pos from phosphopep where id = %s""",(key,))
		x=c.fetchall()[0]
		positions.append(x)
		d[key]=x
	positions.sort()
	for pose in positions:
		for k in keys:
			if d[k]==pose:
				new.append(k)
	return new
		
		
		
		
class Protein:
	def __init__(self,pid,PLENGTH,dictionary,phosphodict, expList,othersList,c,canvas ="Canvas",colors=colors,totallen=900):
		self.colors=colors
		self.PLENGTH=PLENGTH
		self.d = dictionary
		self.p=phosphodict
		self.canvas=canvas
		self.expList=expList
		self.othersList=othersList
		self.totallen=totallen
		self.c=c
		self.pid=pid
	# order phosphorylation sites spatially
	def order(self,keys):
		sites = [PSite(self.PLENGTH,self.p[key],self.totallen) for key in keys]
		orders = [site.d.get("site_pos") for site in sites]
		orders.sort()
		new = []
		for item in orders:
			for s in sites:
				if site.d.get("site_pos")==item: new.append(s.d.get("id"))
		return new

	def drawInJS(self,exp_id):
		offset = 75
		yTop=0
		yBottom = 200+offset
		print '<script type="text/javascript">'
		#print "<!--"
		print 'var jg = new jsGraphics("%s");' % self.canvas
		print 'jg.fillRect(0,'+str(offset+80)+',' +str(self.totallen)+ ',4);'

		# code for spacing names
		if len(self.expList)>0:
			topInterval = 900/(len(self.expList))
		else:
			topInterval = 0
		if len(self.othersList)>0:
			bottomInterval = 900/(len(self.othersList))
		else:
			bottomInterval=0
		

		# end code for spacing names
		# make list of ordered keys for each type of phosphosite
		phosphosites=order(self.p.keys(),self.c)
		pTop=[]
		for item in phosphosites:
			if item in self.expList: pTop.append(item)
#		pTop = [item for item in phosphosites if item in self.expList]
		pBottom = [item for item in phosphosites if item in self.othersList]
		# for each phosphosite in others:
		for phosphosite in phosphosites:
			if phosphosite in self.expList: origin = 'exp'
			else: origin = 'all'
			phos = PSite(self.PLENGTH,self.p[phosphosite],self.totallen,origin =origin)
			print 'jg.setColor("' + phos.getColor() + '");'
			print 'jg.setStroke(3);'
			rcoords=phos.getLineCoords()
			u='%sbrowse/protein.cgi?exp_id=%s&protein_id=%s#' % (urlpath,exp_id,str(self.pid))
			maxnum = 1
			if phos.origin == 'exp':
				if len(self.expList)>maxnum:
					index= pTop.index(phosphosite)
					print "jg.drawLine("+str(rcoords[0])+','+str(rcoords[1]+offset)+','+str(index*topInterval)+','+str(yTop+3+offset)+");"
					print "jg.drawEllipse("+str(index*topInterval)+','+str(yTop+offset)+','+'4,4);'
					print 'jg.setStroke(3);'
					print 'jg.drawLine(' + str(rcoords[0])+','+str(rcoords[1]+offset)+','+str(rcoords[2])+','+str(rcoords[3]+offset)+');'
					t=phos.getText()
					text = ''
					for item in t:
						text = text + item + '<br>'
					t=u+t
                                        print """jg.drawString("<a href='"""+str(t)+"""' style='text-decoration:none'>"""+ text +'<\/a>",'+ str(index*topInterval)+ ',0);' 
				else:
					print 'jg.drawLine(' + str(rcoords[0])+','+str(rcoords[1]+offset-35)+','+str(rcoords[2])+','+str(rcoords[3]+offset)+');'
					ecoords = phos.getEllipseCoords()
					print 'jg.drawEllipse(' + str(rcoords[0]-2)+','+str(offset+ecoords[1])+',4,4);'
					print 'jg.setColor("black");'
					t=phos.getText()
					text = ''
					for item in t:
						text = text + item + '<br>'
					t=u+t
					print """jg.drawString("<a href='"""+str(t)+"""' style='text-decoration:none'>"""+text +'<\/a>",'+str(rcoords[0])+',0);'
																	
														
			else:
				if len(self.othersList)>maxnum:
					index=pBottom.index(phosphosite)
					print "jg.drawLine("+str(rcoords[2])+','+str(rcoords[3]+offset)+','+str(index*bottomInterval)+','+str(yBottom+3)+");"
					print "jg.drawEllipse("+str(index*bottomInterval)+','+str(yBottom)+','+'4,4);'
					print 'jg.drawLine(' + str(rcoords[0])+','+str(rcoords[1]+offset)+','+str(rcoords[2])+','+str(rcoords[3]+offset)+');'
					t=phos.getText()
                                        text = ''
                                        for item in t:
					       text = text + item + '<br>'
					t=u+t
                                        print """jg.drawString("<a href='"""+str(t)+"""' style='text-decoration:none'>"""+text +'<\/a>",'+str(index*bottomInterval)+','+ str(yBottom+5)+');'
					
				else:
					print 'jg.drawLine(' + str(rcoords[0])+','+str(rcoords[1]+offset)+','+str(rcoords[2])+','+str(35+rcoords[3]+offset)+');'
					ecoords = phos.getEllipseCoords()
					print 'jg.drawEllipse(' + str(rcoords[0]-2)+','+str(offset+ecoords[1])+',4,4);'
					t=phos.getText()
					text = ''
					for item in t:
						text = text + item + '<br>'
					t=u+t
					print """jg.drawString("<a href='"""+str(t)+"""' style='text-decoration:none'>"""+text +'<\/a>",'+str(rcoords[0])+','+str(yBottom-20)+');' 
			
		for domain in self.d.keys():
			#print 'jg.setColor("green");'
			dom = Domain(self.d[domain],self.PLENGTH,self.totallen,40)
			coords = dom.getTextRect()
			name = coords[0]
			text= '('+str(dom.d.get("start"))+'-'+str(dom.d.get("stop"))+')'
			#print 'jg.setColor("black");'
			#print 'jg.drawStringRect("'+coords[0]+'",'+str(coords[1])+','+str(40+offset)+','+str(coords[3])+',"'+str(coords[4])+'");'
			#print 'jg.setColor("green");'
			print 'jg.setColor("%s");'%self.colors.get(name,'red')
			print 'jg.fillRect('+str(coords[1])+','+str(60+offset)+','+str(coords[3])+','+'40);'
			print 'jg.setColor("black");'
			print 'jg.drawStringRect("'+text+'",'+str(coords[1])+','+str(70+offset+10+40)+','+str(max(coords[3],len(text)*7))+','+'40);'
			if "Pfam-B" in coords[0]:
			    print 'jg.drawStringRect("'+ "<a target='_blank' href='http://pfam.sanger.ac.uk//pfamb/" + coords[0]+ "'><b>" +coords[0]+'<\/b><\/a>'+'",'+str(coords[1])+','+str(40+offset+23+40)+','+str(max(coords[3],len(coords[0])*8))+',"'+str(coords[4])+'");'
			else:
			    print 'jg.drawStringRect("'+ "<a target='_blank' href='http://pfam.sanger.ac.uk/family?entry=" + coords[0]+ "'><b>" +coords[0]+'<\/b><\/a>'+'",'+str(coords[1])+','+str(40+offset+23+40)+','+str(max(coords[3],len(coords[0])*8))+',"'+str(coords[4])+'");'

			
		print 'jg.paint();'
		#print "-->"
		print '</script>'
		
# origin is "exp" if it comes from the particular experiment of interest
class PSite:
	def __init__(self,PLENGTH,dictionary,BAR,origin="exp"):
		self.PLENGTH=PLENGTH
		self.d=dictionary
		self.BAR=BAR
		self.origin=origin
	
	def convert(self,xvalue):
		return int(1.0*xvalue/self.PLENGTH*self.BAR)
	def getLineCoords(self):
		x1=x2=self.convert(self.d.get("site_pos"))
		if self.origin =='exp':
			off=30
			y1=0+off+10
			y2=50+off
		else:
			off=30
			y1=50+off
			y2=100+off
		return x1,y1,x2,y2
	def getText(self):
		return self.d.get("name")
	def getColor(self):
		value = self.d.get("name")[0]
		colors = {"S":"ccffcc","Y":"ffcc66","T": "cccc66","K": "cfcfcf"}
		
		return colors.get(value,"ccffcc")
	def getEllipseCoords(self):
		if self.origin == 'exp': return (self.getLineCoords()[0],0,4,4)
		else: return (self.getLineCoords()[0],160,4,4)


