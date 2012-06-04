#!/usr/local/bin/python
from pathname import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import sys

sys.path.append(path)
sys.path.append(path+'includes')
from template import *
from smtplib import SMTP,SMTPSenderRefused
def main():
	form = cgi.FieldStorage()
	expid = form.getvalue("expid",None)
	db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
	c=db.cursor()
	printHeader("Summary: Error Log", expid)
	printExpHeader(expid, c);
	logtext = getErrorLog(expid,c)
	printSummary(logtext)
	datatext = getData(expid,c)
	if logtext=='': print "<BR><BR>Invalid experiment id<BR><BR>"
	else: printErrors(logtext,datatext)
	printFooter()

def printSummary(logtext):
	accepted = logtext.count("ACCEPTED")
	rejected = logtext.count("REJECTED")
	print "<BR><BR><table>"
	print "<tr><td>Data load summary:</td><td></td></tr>"
	print "<tr><td>Total Peptides: </td><td>%s</td></tr>"%(accepted+rejected)
	print "<tr><td>Accepted: </td><td>%s</td></tr>"%accepted
	print "<tr><td>Rejected:</td><td>%s</td></tr>"%rejected
	print "</table>"
	print "<BR><BR>"
def concat(list,delimiter):
	string = ''
	for item in list: string+=item+delimiter
	return string
def getData(expid,c):
	if expid ==None: return ''
	query = """select dataset from experiment where id = %s"""%expid
	c.execute(query)
	x=c.fetchall()
	if len(x)>0:
		logfile = x[0][0].strip()
		try:
			f = open(dataPath+logfile,'r')
			text = f.read()
			### check for obscure error with blank lines?
			text = text.replace('\n\n','\n')
			f.close()
			return text
		except IOError:
			print "Error reading data file for accessions."
			return ''
				

	return ''

def getErrorLog(expid,c):
	if expid ==None: return ''
	query = """select errorLog from experiment where id = %s"""%expid
	c.execute(query)
	x=c.fetchall()
	if len(x)>0:
		logfile = logPath+x[0][0].strip()
		try:
			f = open(logfile,'r')
			text = f.read()
			f.close()
			return text
		except IOError: return ''
	return ''

def printErrors(text,datatext):
	data = text.split('--------------------------------------------------------------')
	dataLines = []
	rejected = []
	for item in data:
		if item[0:5].strip()=="Data" and '---' not in item:
			dataLines.append(item)
	for item in dataLines:
		if "ACCEPTED" not in item:
			rejected.append(item)
	if len(rejected)== 0:
		return 
	print "<BR><BR><div style='width:80%;'><table>"
	print "<tr><td><h4>The following errors occurred while loading your dataset:</h4></td></tr>"

	
	for item in rejected:
		if testError(item):
			lines = item.strip().split('\n')
			try:
				dataLine = lines[0]
				print "<tr><td>"+dataLine+"</td></tr>"
			except IndexError: pass
			peptide=getPeptide(datatext,int(dataLine.split('# ')[1].strip()))

			accession = getAcc(datatext,int(dataLine.split('# ')[1].strip()))
			
			lines = sets.Set(lines)
			for line in lines:
				print "<tr><td>"
				if line[0:5]=="ERROR": translateError(line,text,peptide,accession)
			print "</td></tr>"
			print "<tr><td>"+"_"*100+"</td></tr>"
	
	print "</table></div><BR><BR>"

def testError(rejected):
	lines = rejected.strip().split('\n')
	dataLine = lines[0]
	return True
def translateError(error,text,peptide,accession):
	url_dict={'swissprot':"http://ca.expasy.org/uniprot/",'entrez_protein':'http://www.ncbi.nlm.nih.gov/protein/','gi':'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=','GO':'http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term='}
	type=returnAccType(accession.strip())
	if type=='':
		type="undefined"
	if "undefined" in type:
		accession = "<a target='_blank' href='%s'>%s</a>"%(url_dict['entrez_protein']+accession,accession)
		print "<br>Protein could not be loaded. Accession %s is of an unrecognized or unsupported type. See supported accession types on our <a target='_blank' href='%stitle=Main_Page#What_accession_types_are_supported.3F'>help</a> page.<br>" % (accession,helpPath)
	## link accession to entrez
	elif 'gi' in type:
		accession = "<a target='_blank' href='%s'>%s</a>"%(url_dict['gi']+accession,accession)
	elif 'entrez_protein' in type:
		accession = "<a target='_blank' href='%s'>%s</a>"%(url_dict['entrez_protein']+accession,accession)
		
	elif 'swissprot' in type:
		accession = "<a target='_blank' href='%s'>%s</a>"%(url_dict['swissprot']+accession,accession)
	
	else:
		pass
	if 'No pfam predictions exist' in error:
		print "<font color='red'>ERROR: </font>Internal Server error (missing domains)"
		## email admin_db
		#x=SMTP()
		#x.connect("euphrates.mit.edu")
		#try:
			#x.sendmail("ptmscout_admin@mit.edu",'boots1602@gmail.com',"(Kristen, I am trying out the load data error log parser part that sends an email to the db admin)\nPTMScout database load error: No pfam predictions:\n" +error)
		#except SMTPSenderRefused: pass
	
		
	elif 'More than one pfam prediction' in error:
		print "<font color='red'>ERROR: </font>Internal Server error (missing domains)"
		#print "<font color='red'>ERROR: </font>Internal Server error (more than one pfam site domain prediction)"
		#x=SMTP()
		#x.connect("euphrates.mit.edu")
		#try:
			#x.sendmail("ptmscout_admin@mit.edu",'boots1602@gmail.com',"(Kristen, I am trying out the load data error log parser part that sends an email to the db admin)\nPTMScout database load error: More than one pfam prediction:\n" +error)
	
		#except SMTPSenderRefused: pass
		## email admin_db
	elif 'Protein does not exist in database' in error:
		print "Problem loading protein accession: %s"%accession
		acc= accession.split('>')[0].split('<')[0].strip()
		errorLine=''
		for item in text.split('\n'):
			if acc in item:
				errorLine = item
		link = ''
		
		possibleErrors = [line for line in text.split('\n') if accession in line]
		for line in possibleErrors:
			if 'Accession does not lead to protein sequence' in line:
				print "<br>Protein could not be loaded. Accession %s leads to a non-protein entry."%accession
			if "can't match stream return to accession" in errorLine or "object in batch query" in errorLine:
				print "<br>A record for %s could not be retrieved"%accession
		print "<br>"
	elif 'No phosphorylation site' in error or 'cannot find phosphorylation in peptide' in error:
		#print error
		print "<br>Peptide %s does not have a phosphorylation site. Please indicate the phosphorylation site with a lowercase s, t, or y, or preceed it with a p (i.e. pS, pT, pY). Please see our <a target='_blank' href='%stitle=Main_Page#What_accession_types_are_supported.3F'>help</a> page for dataset standards."%(peptide,helpPath)
	elif 'Cannot find peptide in sequence' in error:
		#accession = ''
		print "<br>Peptide %s could not be found for protein loaded in database with accesion %s. Please review peptide for bad terminology (a non-amino acid character) and check that %s is correct.<br>" % (peptide, accession, accession)
	elif "Problem during alignment" in error:
		print "<BR>Problem during alignment of peptide %s. Check accession %s and check that peptide contains only known amino acid codes." %(peptide,accession)
	elif "find peptide in sequence" in error:
		print """Peptide %s could not be found for protein loaded in database with accession %s"""%(peptide, accession)
	else:
		## email me
		print "Unknown error. The error has been sent to the site administrator."
		x=SMTP()
		x.connect("localhost")
		try:
			x.sendmail(sysEmail,sysEmail,"PTMScout database load error: " +error)
		except SMTPSenderRefused: pass

	
def returnAccType(accession):
	command = """perl -I"%s" entrez.pl '%s'"""%(libPath,accession)
	f=os.popen(command)
	return f.read()
def getAcc(text,linenumber):
	text = text.replace('\r','\n')
	text = text.replace('\n\n','\n')
	columns = [item.strip() for item in text.strip().split('\n')[0].split('\t')]
	accCol=None
	for item in columns:
		if re.match("acc", item.lower()):
			accCol = columns.index(item)
	if accCol==None: return''
	else:
		try:
			return text.split('\n')[linenumber-1].split('\t')[accCol]
		except IndexError: return ''
	
def getPeptide(text, linenumber):
	text = text.replace('\r','\n')
	text = text.replace('\n\n','\n')
	columns = [item.strip() for item in text.strip().split('\n')[0].split('\t')]
	pepCol=None
	
	for item in columns:
		if re.match("pep", item.lower()) != None or re.match("pep:tryps",item.lower())!=None:
			pepCol = columns.index(item)
	if pepCol==None: return''
	else:
		try:
			#print linenumber,text.split('\n')[linenumber-1]
			return text.split('\n')[linenumber-1].split('\t')[pepCol]
		except IndexError: return ''
	

main()	
