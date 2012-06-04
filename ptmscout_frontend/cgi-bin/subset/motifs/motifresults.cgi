#!/usr/bin/python
from pathname import *
import random
import cgi
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import cStringIO
import MySQLdb
import os
import cgitb
if displayPythonErrors:
	cgitb.enable()
import sys
sys.path.append(path+'includes')

from sets import Set
import time
import re
from template import *

def main():
	print "<html><head><title>PTMScout - Motif Tool Results</title></head><body>"
#	print "Results loading<br>"
	
	
	db = MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
	c=db.cursor()
	form = cgi.FieldStorage()
	expid = form.getvalue('expid','7')
	printHeader("Motif Results", expid);
	include = form.getvalue('include','off')
	type = form.getvalue('motif_type','peptide')
	logfilename = motifPath+'logfile.txt'
	numfilename = motifPath+'numsfile.txt'
	printExpHeader(expid,c)
#	print "DEBUG: "+motifPath+"<br>"

#	os.system('ls '+motifPath)
#	print "<br>"
	timer = int(time.time())
	out = ''

	## get email address if user wants results emailed
	email = form.getvalue("emailresults","off")
	emailaddress = form.getvalue("emailaddress"," ")
	if invalid(emailaddress): emailaddress = " "
## 	########################################
## 	### get number, make sure it's not in use
	try:
		try:
			f = open(numfilename,'r')
		except IOError:
			os.system('touch '+numfilename)
			os.system('chmod 777'+numfilename)
			f=open(numfilename,'r')
		#print "DEBUG: opening numfilename<br>"
		lines = f.readlines()
		nums = [int(line.split('\t')[0].strip()) for line in lines if line.split('\t')[0].strip().isalnum()]
		num = random.randint(0,100000)
		#print "DEBUG: num before is "+str(num)+"<br>"
		while num in nums:
			num = random.randint(0,100000)
			#print "DEBUG: num in loop is "+str(num)+"<br>"
		f.close()

		## write that number back to the number file so we know it's in use now
		f = open(numfilename,'a')
		f.write(str(num)+'\t'+str(timer)+'\t'+emailaddress+'\n')
		f.close()
	except IOError:
		print "Error generating key"
		num = random.randint(0,100000)
	url = "%ssubset/motifs/motifFetch.py?num=%s"%(urlpath,num)
	#use this random key in the log generation.
	print "Motif results can be found at the following url for the next 24 hours (it may take several minutes for the results to be available at this url):<br><a href='%s'>%s</a><BR><BR>"%(url,url)
	#### rewrite query file with new number
	try:
		querynum = form.getvalue("querynum",None)
		querystr = open('%squery%s'%(outPath,querynum),'r').read()
		open('%squery%s'%(outPath,num),'w').write(querystr)
		print querystr
	except IOError:
		print "<BR><BR>Error retrieving foreground search<BR><BR>"
		querystr = ""
## 	#########################################

	try:
		tablefile = open('%ssubset/motifs/tables/table%s.txt'%(path,num),'w')
	except IOError:
		print "<BR>error openingtablefile<BR>"
		tablefile = open('/tmp/ptmscout/table.txt','w')
	tablefile.write('Motif\tFG\tBG\tP-value\n')
## 	######################################
## 	### get foreground and background
	if type=='peptide':
		foreground = form.getvalue('foreground','')
	else:
		foreground = form.getvalue('fgprotein','')
	
	if foreground !='':
		foreground = foreground.split(':')[:-1]

	if type=='peptide':
		background = form.getvalue('background','')
	else:
		background = form.getvalue('bgprotein','')
	if background !='':
		background = background.split(':')[:-1]
## 	###########################################

## 	#########################################
## 	### write foreground and background to file
	
	### for peptide:
	if type=='peptide':
		fg,bg=[],[]
		for item in foreground:
			query="""select pep_aligned from phosphopep where id = %s"""%(item)
			c.execute(query)
			x=c.fetchall()
			if x!=None:
				fg.append(fifteenmerize(x[0][0]))
		for item in background:
			query="""select pep_aligned from phosphopep where id = %s"""%(item)
			c.execute(query)
			x=c.fetchall()
			if x != None:
				bg.append(fifteenmerize(x[0][0]))
		os.chdir(motifPath)
		fgfilename = motifPath+'foreground%s'%num
		bgfilename= motifPath+'background%s'%num
		try:
			f=open(fgfilename,'wb')
			for item in fg[:-1]:
				f.write(item+ '\n')
			f.write(fg[-1])
			f.close()
		except IOError: print "<h3>Error writing foreground file. If problem continues email site administrator</h3>"
		except IndexError:
			print "<h3>Error. Empty foreground</h3>"
			printFooter()
			sys.exit(0)
		try:
			f=open(bgfilename,'wb')
			for item in bg[:-1]:
				f.write(item+'\n')
			f.write(bg[-1])
			f.close()
		except IOError: print "<h3>Error writing background file. If problem continues email site administrator</h3>"
		numpept = max(2,int(len(fg)*0.05))
		numfg = len(fg)
		numbg = len(bg)
	## for protein:
	else:
		fg,bg={},{} # want to make dict[protein number]=list of peptides
		for item in foreground: # item is a protein
#				query="""select pep_aligned from phosphopep where protein_id = %s"""%(item)
                                if include=='on':
					query="""select pep_aligned from phosphopep where protein_id = %s group by pep_aligned"""%item
				else:
					query = """select pep_aligned from phosphopep join MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id on MS_phosphopep.phosphopep_id = phosphopep.id where experiment_id = %s and MS.protein_id = %s group by pep_aligned"""%(expid,item)
					
				peplist = []
				c.execute(query)
				x=c.fetchall()
				if x != None:
					for pep in x:
						peplist.append(pep[0])
					fg[item]=peplist
		for item in background:
				if include=='on':
					query="""select pep_aligned from phosphopep where protein_id = %s group by pep_aligned"""%(item)
				else:
				        query = """select pep_aligned from phosphopep join MS_phosphopep join MS on MS.id = MS_phosphopep.MS_id on MS_phosphopep.phosphopep_id = phosphopep.id where experiment_id = %s and MS.protein_id = %s group by pep_aligned"""%(expid,item)
				peplist = []
				c.execute(query)
				x=c.fetchall()
				if x != None:
					for pep in x:
						peplist.append(pep[0])
					bg[item]=list(Set(peplist))
		fgfilename = motifPath+'foreground%s'%num
		bgfilename= motifPath+'background%s'%num
		try:
				f=open(fgfilename,'wb')
				for item in fg:
					query = """select acc_gene from protein where id = %s""" % item
					c.execute(query)
					x=c.fetchall()
					if len(x)>0:
						acc = x[0][0]
					else:
						acc = item
					f.write('>'+acc+ '\n')
					for pep in fg[item]:
						f.write(pep+'\n')
			
				f.close()
		except IOError: print "<h3>Error writing foreground file. If problem continues email site administrator</h3>"
		try:
				f=open(bgfilename,'wb')
				for item in bg:
					query = """select acc_gene from protein where id = %s""" % item
					c.execute(query)
					x=c.fetchall()
					if len(x)>0:
						acc = x[0][0]
					else:
						acc = item
					f.write('>'+acc+'\n')
					for pep in bg[item]:
						f.write(pep+'\n')
				f.close()
		except IOError: print "<h3>Error writing background file. If problem continues email site administrator</h3>"
			
		numfg = len(fg)
		numbg = len(bg)


############# if emailing results####################3
	if email == 'on':
		## check if email address is valid
		## print appropriate message
		if not (invalid(emailaddress)):
			print "Results will be sent to %s within the next 24 hours" % emailaddress
		else:
			print "%s is an invalid email address" % emailaddress
		print "<BR><BR>"
		## exit
		
######################################################


## 	#########################################

## 	#########################################
## 	##################submit job to server###
## 	#########################################


## 	#######################################
## 	### calculate results
	
	os.chdir(motifPath)
	if type=="peptide":			
			command1 = "perl c_stat_sig_peptide_refpass_depth_memory.phospho_acetyl.pl --foreground %s --background %s --numpept %s " % (fgfilename.split('/')[-1],bgfilename.split('/')[-1],max(2,numpept))
			

	else:
			command1 = "perl c_stat_sig_protein_refpass_depth_memory.phospho_acetyl.pl --foreground %s --background %s" % (fgfilename,bgfilename)
	os.popen(command1)
	command2="%sprocess_motifs.pl %slog.%s > %sweboutput%s.p" % (motifPath, motifPath,fgfilename.split('/')[-1]+'.'+bgfilename.split('/')[-1],motifPath,num)
	os.popen(command2)

	outputfile = motifPath+'weboutput%s.p'%num
	ready = True
	try:
		f=open(outputfile,'rb')
		output = f.read()
		f.close()
	except:
		ready = False
		#here add call to create results
		print "<BR>Results not ready yet. Check the url above, and results will be posted within the next hour.<BR>"
		output = ''
	output = output.split('\n')
	data = []
	# sort output on p-value
	for item in output:
		if len(item)>1:
			data.append(  (float(item.split('|')[4].strip()),item))
	data.sort()
## 	########################################



## 	############################################
## 	## display results
	out+= '<a href="%ssubset/motifs/makeforeground.txt?num=%s">Foreground</a>'%(urlpath,num)+'<br>'
	out+= '<a href="%ssubset/motifs/makebackground.txt?num=%s">Background</a>'%(urlpath,num)+'<br>'
	if len(data)>0:
		out += '<center><a href="%ssubset/motifs/table.txt?num=%s">Export results to file</a></center>'%(urlpath,num)
		out+= '<center><table   width = 60%% >'
		if type == "peptide":
			out+= '<tr bgcolor="cccccc"><th>Motif</th><th>In Foreground</th><th>In Background</th><th>P-value</th></tr>'
		else:
			out+= '<tr bgcolor="cccccc"><th>Motif</th><th>In Foreground</th><th>In Background</th><th>P-value</th></tr>'
		for item in data:
			formattedP = "%.2e" % item[0]
			base, power = formattedP.split('e')
			pval = base+' x ' + '10<sup>%s</sup>' % str(int(power))
			realmotif = item[1].split('|')[1]
			motif = realmotif
			inFG = item[1].split('|')[2].split('/')[0].strip()
			inBG = item[1].split('|')[3].split('/')[0].strip()
			out+= '<tr align="center">'
			out+= '<td><tt>%s</tt></td><td>%s</td><td>%s</td><td>%s</td>' %(realmotif.strip(),inFG+'/'+str(numfg),str(inBG)+'/'+str(numbg),pval)
			tablefile.write(realmotif.strip()+'\t'+str(inFG)+'/'+str(numfg)+'\t'+str(inBG)+'/'+str(numbg)+'\t'+str(item[0])+'\n')
			out+= '</tr>'
		out+= '</table></center>'
	else:
		if ready: print "<h3>No results yet.</h3>"
	out+= "<BR>"*5		
	fgf = fgfilename.split('/')[-1]
	bgf = bgfilename.split('/')[-1]
	intermediate = 'weboutput%s.p.s'%num
	intermediate2 = 'log.%s.%s'%(fgf,bgf)
	command4 = "rm %s" %fgf
	command5 = "rm %s" %bgf
	command6 = "rm %s" %outputfile
	#command7 = "rm %s" %(fgf+'.'+bgf)
	command8 = "rm %s" %intermediate
	command9 = "rm %s" %intermediate2
	#os.popen(command4)
	#os.popen(command5)
	#os.popen(command6)
	os.popen(command8)
	os.popen(command9)

	tablefile.close()
	
	out+= "</body></html>"
	print out
	out =  "<html><head><title>PTMScout - Motif Tool</title></head><body>"+out
	printFooter()

	########################################
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
		middle = peptide.find(residue)
		if middle == 7:
			if len(peptide) <15:
				return peptide+(15-len(peptide))*' '
			else:
				return peptide[0:15]
		elif middle < 7:
			pep= ' '*(7-middle)+peptide
			if len(pep)>15: return pep[0:15]
			else: return pep + (15-len(pep))*' '
		else:
			pep = peptide[middle-7:]
			if len(pep)>15: return pep[0:15]
			else:
				return pep + (15-len(pep))*' '
			
		
	return peptide

GENERIC_DOMAINS = "aero", "asia", "biz", "cat", "com", "coop", \
        "edu", "gov", "info", "int", "jobs", "mil", "mobi", "museum", \
        "name", "net", "org", "pro", "tel", "travel"

def invalid(emailaddress, domains = GENERIC_DOMAINS):
        """Checks for a syntactically invalid email address."""

        # Email address must be 7 characters in total.
        if len(emailaddress) < 7:
            return True # Address too short.

        # Split up email address into parts.
        try:
            localpart, domainname = emailaddress.rsplit('@', 1)
            host, toplevel = domainname.rsplit('.', 1)
        except ValueError:
            return True # Address does not have enough parts.

        # Check for Country code or Generic Domain.
        if len(toplevel) != 2 and toplevel not in domains:
            return True # Not a domain name.

        for i in '-_.%+.':
            localpart = localpart.replace(i, "")
        for i in '-_.':
            host = host.replace(i, "")

        if localpart.isalnum() and host.isalnum():
            return False # Email address is fine.
        else:
            return True # Email address has funny characters.
main()



### on load, motifresults does its normal job
### saves [num] to numfile with timestamp
### run cron job that reads num file, makes outputHTMLS, deletes all numbers older than 24 hours


### to make:
### motifCleanup.py runs every minute:
   ## reads num file
   ## deletes numbers older than 24 hours (also deletes outputHTML,foreground and background files. all else should be deleted by motifresults.cgi)
   ## for remaining numbers: checks if there is an outputHTML for it
   ## if not, creates it
