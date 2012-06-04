#!/usr/local/bin/python
from pathname import *
import cgi
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
import MySQLdb
import sys

sys.path.append(path)
sys.path.append(path+'includes')
import cgitb
if displayPythonErrors:
	cgitb.enable()
from template import *
import copy
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import random
import crypt
from MimeWriter import MimeWriter
import re
import os

try:
	from cStringIO import StringIO
except:
	from StringIO import StringIO
from smtplib import SMTP,SMTPSenderRefused
import urllib
def main():
	db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
	c=db.cursor()
	returnEmail = sysEmail#"ptmscout_admin@mit.edu"
	form = cgi.FieldStorage()
	formdata= [(item,form.getvalue(item)) for item in form.keys() if item != "dataFile"]
	formdata = cgi.urllib.urlencode(formdata)
	printHeader("PTMScout - Load Dataset")
	query = [(item, form.getvalue(item,'')) for item in form.keys() if item !="dataFile"]
	querystring= urllib.urlencode(query)
	formdata = ''
	for item in form.keys():
		if item !="dataFile":
			key = item
			val = form.getvalue(item,'')
			formdata+="<input type=hidden name='%s' value='%s'>"%(key,val)
	if (not(form.has_key('dataFile')) or form.getvalue('dataFile','')=='') and not(form.has_key('check')):
		print """<br><br>No input file loaded."""
		

	else:
		### no matter what load a copy to tmp
		file = form.getvalue('dataFile','')
		goodToGo=True
		if not checkForm(form,c):
			goodToGo = False
		if form.has_key('dataFile'):
			tmpfilename = '/tmp/'+form['dataFile'].filename+str(random.randint(0,10000000))
		else:
			tmpfilename = '/tmp/temp.txt.try'+str(random.randint(0,10000000))
		try:
			## write temp file with empty lines removed
			open(tmpfilename,'wb').write(removeEmptyLines(file))
			## convert newlines to UNIX
			os.system("dos2unix -n %s %s"%(tmpfilename,tmpfilename))
			os.system("mac2unix -n %s %s"%(tmpfilename,tmpfilename))
		except IOError:
			print "<BR>Error writing datafile. Contact administrator.<BR>"
			goodToGo = False
	        ### load their file to special directory
		if not(form.has_key('check')):
			file = form.getvalue('dataFile','')
			filename = dataPath+form['dataFile'].filename+str(random.randint(0,10000000))
			try:
				open(filename,'wb').write(file)
			except IOError:
				print "<br><br>Error loading file. Contact administrator<br><br>"
				goodToGo = False
		else:
			try:
				filename = form.getvalue("filename","")
				file=open(filename,'r').read()
			except IOError:
				#print filename
				print "<BR><BR>Error reading data file. Contact administrator<Br><BR>"
				goodToGo = False
		name = form.getvalue("name","")
		email = form.getvalue("youremail","")
		expemail = form.getvalue("email","")
		expname = form.getvalue("expname","")
		authors = form.getvalue("authors","")
		description = form.getvalue("description","")
		pmid = form.getvalue("pmid",0)
		loadpm = form.getvalue("loadpm","off")
		extendChange = form.getvalue("extendChange","")
		if extendChange != "": expname = "[%s]%s"%(extendChange,expname)
		if pmid == "" or pmid == " ": pmid = 0
		url = form.getvalue("url","NA")
		if url == "" or url == " ": url = "NA"
		if url == "NA" and pmid != 0:
			url = "http://www.ncbi.nlm.nih.gov/sites/entrez?holding=&db=pubmed&cmd=search&term="+pmid
		published = form.getvalue("published","NO")
		journal = form.getvalue("journal","")
		volume = form.getvalue("volume","")
		pages = form.getvalue("pages","")
		date = form.getvalue("date","")
		ambiguity = form.getvalue("ambiguity","NO")		
		notes = form.getvalue("notes","")
		export = form.getvalue("export","YES")
		loadType = form.getvalue("loadType","new")
		if loadType == "reload":
			replace = form.getvalue("change","")
			if replace == "": loadTYpe = "new"
		else: replace = ""
		if loadType == "exportLink":
			link = form.getvalue("link","")
			link = getRealParent(link,c)
		else: link = "0"
			
		### get modification type
		mods = []
		if form.has_key("modT"): mods.append("T")
		if form.has_key("modY"): mods.append("Y")
		if form.has_key("modS"): mods.append("S")
		if form.has_key("modK"): mods.append("K")
		modString = ''
		for mod in mods:
			modString += mod+","
		try:
			if modString[-1]==',': modString= modString[:-1]
		except IndexError:
			print "Error: No modification type designated."
			goodToGo = False
		
		### preliminary tests on file
		checkData = form.getvalue('check',True)
		if checkData == "False": checkData = False
		if checkData == "True": checkData = True

		argDict = {"YES":"1", "NO": "0"}
		if goodToGo:
			if prelimTestFile(file,tmpfilename,formdata,checkData = checkData, query=querystring,filename=filename):
				export = "NO"
				experiment_id_link = "None"
				text = """PTMScout data load\n\n"""
				text+= """User Information:\n"""
				text += """Name: %s\n"""%name
				text += """Email: %s\n\n"""%email
				if loadType == "new" or loadType == "exportLink":
					text += "New dataset:\n\n"
				else:
					text += "Delete experiment %s\n"%replace
				text += filename + ' ' + "'"+expname+"'"  + ' ' + "'"+authors+"'" + ' '+"'"+description+"'"+' '+"'"+expemail+"'" + ' '+str(pmid) + ' '+"'"+ url +"'"+ ' '+argDict.get(published.upper(),'0') + ' '+ argDict.get(ambiguity.upper(),'0') + ' '+"1" + ' '+ link + ' '+"'" + email+"'"+' ' + "'" + modString+"'" + " '"+journal+"'" + " '" + date+ "'" + " '" + str(volume)+"'"+ " '"+pages + "'" + " '"+PERL_DB+"'"
				
				text+='\n\n'


				text += """Notes: %s\n\n"""%notes
				tempfile = StringIO()
				mw = MimeWriter(tempfile)
				mw.addheader('to',returnEmail)
				mw.addheader('from',email)
				mw.addheader('subject','[PTMScout] Load dataset')
				mw.startmultipartbody('mixed')
				sw=mw.nextpart()
				j = sw.startbody('text/plain')
				j.write(text)
				mw.lastpart()
				message = tempfile.getvalue()
				mailer = SMTP()
				mailer.connect("localhost")

				try:
					mailer.sendmail(email,returnEmail,message)
					print "<br><br>Dataset submitted successfully. You will receive an email notification when your file has been processed<br><br>"
				except SMTPSenderRefused:
					print "<br><br>Invalid email address. Please correct and try again.<br><br>"
				mailer.close()
				

			
		else:
			print """<BR><BR>File not loaded. Error reading data file<BR><BR><BR>"""

	printFooter()



def prelimTestFile(file,tmpfilename,formdata,checkData = True,query='',filename=''):
	"""prelimTestFile(file,tmpfilename,formdata,checkData=True,query='',filename='')
	returns True if the file is acceptable, False and an error message if not"""
	file = removeEmptyLines(file)
	lines = file.split('\n')

	columns = [item.strip().replace("\"","") for item in lines[0].split('\t')]
	### check for accession column and peptide column
	accCol = None
	pepCol = None
	dataCols=[]
	runCol = None
	for item in columns:
		if re.match("acc",item.lower()) != None or re.match("acc:",item.lower()) != None:
			if accCol != None:
				print """<BR><BR>ERROR: More than one accession column found<BR><BR>"""
				return False ##means we found more than one acc column
			else:
				accCol = columns.index(item)
		if re.match("pep", item.lower()) != None or re.match("pep:tryps",item.lower())!=None:
			if pepCol != None:
				print """<BR><BR>ERROR: More than one peptide column found<BR><BR>"""
			else:
				pepCol = columns.index(item)
		if re.match("data",item.lower()) != None or re.match("data:",item.lower())!=None:
			dataCols.append(columns.index(item))
		if re.match("run",item.lower())!=None:
			runCol = columns.index(item)
	if accCol == None:
		print """<BR><BR>ERROR: No accession column found<BR><BR>"""
		return False ## we found no acc column
	if pepCol == None:
		print """<BR><BR>ERROR: No peptide column found<BR><BR>"""
		return False ## we found no pep column
	accs=[]
	types=[]
	peps = []
	pepFormats = []
	needRun = False
	accPeps = []
	for i in range(1,len(lines)-1):
		try:
			acc = lines[i].split('\t')[accCol].strip().replace('\"','')
			
			if acc == None or acc=="null" or acc=="":
				print """<BR><BR>ERROR (Line %s): Null accession found.<BR><BR>"""%i
				return False
			pep = lines[i].split('\t')[pepCol].strip().replace('\"','')
			if pep == "" or (pep.isalnum() and not(pep.isalpha())):
				print """<BR><BR>ERROR (Line %s): Missing or corrupted peptide<BR><BR>"""%i
				return False
			
			
			if (acc,pep) in accPeps:
				#print (acc,pep)
				needRun = True
			accPeps.append((acc,pep))
			type = returnAccType(acc)
			accs.append(acc)
			types.append(type)
			
			
			peps.append(pep)
			if 'y' not in pep and 't' not in pep and 's' not in pep and 'pY' not in pep and 'pT' not in pep and 'pS' not in pep:
				pepFormats.append(False)
			else:
				pepRight = True
				for letter in pep:
					if letter not in ['y','s','t','p','G','A','V','L','I','M','F','W','P','S','T','C','Y','N','Q','D','E','K','R','H']:
						pepRight = False
				pepFormats.append(pepRight)
			
		except IndexError:
			print """<BR><BR>ERROR (Line %s): Line incorrectly formatted.<BR><BR>"""%i
			return False
	
	##if len(["here" for item in types if 'undefined' in item])>0.5*len(types): # this checked if more than half of accessions bad
	## now sample 20 types, if more than 10 are bad give error
	copyTypes = copy.deepcopy(types)
	if len(copyTypes) > 20: copyTypes = copyTypes[0:20]
	if len(["here" for item in copyTypes if 'undefined' in item])>10:
	
		print "<BR><BR>ERROR: More than half of tested accessions have unknown type.<br><BR>"
		return False
	if pepFormats.count(False)>0.5*len(pepFormats):
		print "<BR><BR>ERROR: More than half of peptides are incorrectly formatted. See documentation for information on correctly formatting files.<BR><BR>"""
		return False

	### check that dataCols formatted correctly
	#print dataCols
	for col in dataCols:
		if not (lines[0].split('\t')[col].count(':')==2):
			print """<BR><BR>ERROR: Data column headers formatted incorrectly. Should be data:type:label<BR><BR>"""
			return False
	if checkData:
		if len(dataCols)==0:
			print """<BR><BR>WARNING: No data columns found.<BR><BR>"""
			print """<BR><BR><form action="load.cgi" method="POST"><input type=hidden value=False name='check'><input type=hidden value=%s name='filename'><input type=submit value="Continue">%s</form><BR><BR>"""%(filename,formdata)
			return False
		else:
			print """<BR><BR>You have entered data columns. <a target = '_blank' href='preview.cgi?file=%s'>Preview data</a> %s"""%(tmpfilename,helper("Load_Dataset#Preview_Data"))
			print """<BR><BR><form action="load.cgi" method="POST"><input type=hidden value=False name='check'><input type=hidden value=%s name='filename'><input type=submit value="Continue">%s</form><BR><BR>"""%(filename,formdata)
			return False
	### check if need run column
	#if needRun and runCol == None:
	if needRun and runCol == None:
		print """<BR><BR>ERROR: You must enter a run column<BR><BR>"""
		print """<table><tr><td colspan='4'>Example:</td></tr>
		<tr><td>Acc</td><td>pep</td><td>run</td><td>data</td></tr>
		<tr><td>Acc1</td><td>pep1</td><td>1</td><td>data</td></tr>
		<tr><td>Acc1</td><td>pep1</td><td>2</td><td>data</td></tr>
		<tr><td>Acc2</td><td>pep1</td><td>1</td><td>data</td></tr>
		<tr><td>Acc3</td><td>pep1</td><td>1</td><td>data</td></tr>
		</table><BR><BR>"""
		return False
	return True

#from perlfunc import perlfunc,perlreq,perl5lib
#@perlfunc
#@perlreq("entrezTools.pm")
def returnAccType(accession):
	#pass
	f=os.popen("""perl -I"%s" entrez.pl '%s'"""%(libPath,accession))
	return f.read()

def getRealParent(expid,c):
	query = """select experiment_id from experiment where id = %s"""%expid
	c.execute(query)
	try:
		x=c.fetchall()[0][0]
		if int(x)==0: return str(expid)
		else: return str(x)
	except IndexError: return str(expid)

def removeEmptyLines(fileText):
	""" removeEmptyLines(fileText) removes empty lines from the text"""
	text = ""
	fileText.replace('\r','\n')
	lines = fileText.split('\n')
	for line in lines:
		if line.strip()!="":
			text+= line.strip()+'\n'
	return text

def getPublicationInfo(form,c):
	""" get publication information from pmid fetch"""
	outputfile = outPath+"pmid"+str(random.randint(0,10000000))
	pmid= form.getvalue("pmid",0)
	os.system("""perl -I"%s" %scall_printRefFromPMID.pl %s %s"""%(libPath,scriptPath,pmid,outputfile))
	try:
		f = open(outputfile,'r')
		text = f.read()
	except IOError:
		return ""
	if text == "": return ""
	journal,date,volume,pages = "","","",""
	try:
		lines = text.split('\n')
		volume = lines[1].split("\t")[1]
		date = lines[2].split("\t")[1]
		pages = lines[3].split("\t")[1]
		journal = lines[5].split("\t")[1]
	except IndexError:
		return ""
	## for line in lines:
## 		item, value = line.split("\t")[0],line.split("\t")[1]
## 		if item.strip() == "journal":
## 			journal = value
## 		if item.strip() == "pub_date":		
## 			date = value
## 		if item.strip() == "pages":
## 			pages = value
## 		if item.strip() == "volume":
## 			volume = value
	#print journal,date,volume,pages
	return journal,date,volume,pages
	
		
		

def checkForm(form,c):
	"""perform checks again on the form, shouldn't rely on javascript to check"""
	###check that required fields are not empty
	if form.getvalue("email","") == "":
		print "<BR><BR>Missing email address<BR><BR>"
		return False
	if form.getvalue("youremail","")=="":
		print "<BR><BR>Missing user email address<BR><BR>"
		return False
	#if form.getvalue("dataFile","")=="":
		#print "<BR><BR>No data file<BR><BR>"
		#return False
	if form.getvalue("expname","")=="":
		print "<BR><BR>Missing experiment name<BR><BR>"
		return False
	if form.getvalue("authors","")=="":
		print "<BR><BR>Missing authors<BR><BR>"
		return False
	if form.getvalue("description","")=="":
		print "<BR><BR>Missing description<BR><BR>"
		return False
	###check that email addresses are formatted correctly
	if not isAddressValid(form.getvalue("email","")):
		print "Invalid email address"
		return False
	if not isAddressValid(form.getvalue("youremail","")):
		print "Invalid email address"
		return False

	###validate agreement checked
	if not form.getvalue("agree","off") == "on":
		print "You must agree to the terms of use."
		return False

	###validate have change
	if form.getvalue("loadType","new") == "exportLink":
		if form.getvalue("extendChange","")=="":
			print "You must enter a description of the change to the dataset."
			return False
	
	###
	## check pmid and publication stuff
	published = form.getvalue("published","NO")
	try:
		pmid = int(form.getvalue("pmid",0))
	except ValueError:pmid = 0
	loadpm = form.getvalue("loadpm","off")
	journal = form.getvalue("journal","")
	volume = form.getvalue("volume","")
	pages = form.getvalue("pages","")
	date = form.getvalue("date","")
	## check that pmid isn't already in database
	query = """select PMID from experiment group by PMID"""
	c.execute(query)
	x=c.fetchall()
	pmids = [int(item[0]) for item in x]
	if pmid in pmids and pmid !=0 and form.getvalue("loadType","new") != "reload" and form.getvalue("loadType","new") == "new":
		print "Experiment already exists, please return to data load and choose to load Extension/Change of an existing dataset if you would like to continue loading a modified version of this experiment."
		printFooter()
		sys.exit(0)
	###
	if published == "YES":
		if pmid == 0 and loadpm == "on":
			print "Invalid pub med ID"
			return False
		elif pmid !=0 and loadpm == "on":
			pass
		elif journal == "" or volume == "" or pages == "" or date == "":
			print "Must either enter publication data or load from pub med ID"
			return False
		
	######
	return True
from string import *

rfc822_specials = '()<>@,;:\\"[]'

def isAddressValid(addr):
    # First we validate the name portion (name@domain)
    c = 0
    while c < len(addr):
        if addr[c] == '"' and (not c or addr[c - 1] == '.' or addr[c - 1] == '"'):
            c = c + 1
            while c < len(addr):
                if addr[c] == '"': break
                if addr[c] == '\\' and addr[c + 1] == ' ':
                    c = c + 2
                    continue
                if ord(addr[c]) < 32 or ord(addr[c]) >= 127: return 0
                c = c + 1
            else: return 0
            if addr[c] == '@': break
            if addr[c] != '.': return 0
            c = c + 1
            continue
        if addr[c] == '@': break
        if ord(addr[c]) <= 32 or ord(addr[c]) >= 127: return 0
        if addr[c] in rfc822_specials: return 0
        c = c + 1
    if not c or addr[c - 1] == '.': return 0

    # Next we validate the domain portion (name@domain)
    domain = c = c + 1
    if domain >= len(addr): return 0
    count = 0
    while c < len(addr):
        if addr[c] == '.':
            if c == domain or addr[c - 1] == '.': return 0
            count = count + 1
        if ord(addr[c]) <= 32 or ord(addr[c]) >= 127: return 0
        if addr[c] in rfc822_specials: return 0
        c = c + 1

    return count >= 1
main()
