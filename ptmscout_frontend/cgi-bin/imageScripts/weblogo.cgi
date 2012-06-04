#!/usr/local/bin/python
from pathname import *
import cgi
import cgitb
if displayPythonErrors:
	cgitb.enable()
import sys,os

sys.path.append(path)
sys.path.append(path+'includes')
import MySQLdb
db = MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
c=db.cursor()
import random
def main():
	form = cgi.FieldStorage()
	filename = form.getvalue("filename",'')
	## expid = form.getvalue("expid",8)
## 	selected = form.getvalue("selectedids","")
## 	filename = "%sexp%s%s.fasta"%(fastaPath,expid,random.randint(0,100000000))
## 	if expid == None or expid == "None":
## 		selected = selected.split(':')
## 		ids = []
## 		if selected !=['']:
## 			for item in selected:
## 				ids.append(item)
## 		msids=ids
## 		makeFasta2(c,msids,filename)
## 	else:
## 		makeFasta(c,expid,filename)
	command = """%s/seqlogo -c -S -F PNG -f %s > %s.png"""%(weblogoPath,filename,filename)
	os.system(command)
	
	f = open(filename+".png",'r')
	f.seek(0)
	data = f.read()
	print "Content-Type: image/png\nContent-Length: %d\n" % len(data)
	print data
	f.close()

## def makeFasta2(c,msids,filename):
## 	f = open(filename,'w')
## 	id = 1
## 	for item in msids:
## 		query = """select pep_aligned from MS join MS_phosphopep on MS.id = MS_phosphopep.MS_id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id = %s"""%item
## 		c.execute(query)
## 		x=c.fetchall()
## 		x = [thing[0] for thing in x]
## 		for pep in x:
## 			if len(fifteenmerize(pep.strip()))==15:
## 				f.write(">%s\n"%id)
## 				f.write(fifteenmerize(pep.strip())+"\n")
## 		id += 1
## 	f.close()
## def fifteenmerize(peptide):
## 	if len(peptide)==15: return peptide
## 	if 'y' in peptide:
## 		residue = 'y'
## 	elif 't' in peptide:
## 		residue = 't'
## 	elif 's' in peptide:
## 		residue = 's'
## 	elif 'k' in peptide:
## 		residue = 'k'
## 	else:
## 		residue = ''
## 	if residue != '':
## 		peptide = peptide.strip()
## 		middle = peptide.find(residue)
## 		if middle == 7:
## 			if len(peptide) <15:
## 				return peptide+(15-len(peptide))*'.'
## 			else:
## 				return peptide[0:15]
## 		elif middle < 7:
## 			pep= '.'*(7-middle)+peptide
## 			if len(pep)>15: return pep[0:15]
## 			else: return pep + (15-len(pep))*'.'
## 		else:
## 			pep = peptide[middle-7:]
## 			if len(pep)>15: return pep[0:15]
## 			else:
## 				return pep + (15-len(pep))*'.'
			
		
## 	return peptide
## def getMSids(c,expid):
## 	query = """select id from MS where experiment_id = %s"""%expid
## 	c.execute(query)
## 	x=c.fetchall()
## 	return [item[0] for item in x]

## def makeFasta(c,expid,filename):
## 	f = open(filename,'w')
## 	query = """select pep_aligned from MS join MS_phosphopep on MS.id = MS_phosphopep.MS_id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where experiment_id=%s"""%expid
## 	c.execute(query)
## 	x=c.fetchall()
## 	peps = [item[0] for item in x]
## 	id = 1
## 	for pep in peps:
## 		if len(fifteenmerize(pep.strip())) == 15:
## 			f.write(">%s\n"%id)
## 			f.write(fifteenmerize(pep.strip())+"\n")
## 			id += 1
## 	f.close()
		
main()
