import os
import sys
from pathname import *
import random
import time
import re

def main():
        out = ''
        num = sys.argv[1]
        fgfilename = motifPath+'foreground'+str(num)
        bgfilename = motifPath + 'background'+str(num)
        numpept = 0
	### get foreground file and background file

        type="peptide"
        ### detect if we need peptide or protein version
        try:
            if '>' in open(bgfilename,'r').read(): type = "protein"
        except IOError: pass
        
	timer = int(time.time())


	## get numfg
        try:
            numfg = len(open(fgfilename,'r').readlines())
            numbg = len(open(bgfilename,'r').readlines())
        except IOError:
            numfg = 0
            numbg = 0



## 	#######################################
## 	### calculate results
	os.chdir(motifPath)
	if type=="peptide":			
			command1 = "perl c_stat_sig_peptide_refpass_depth_memory.phospho_acetyl.pl --foreground %s --background %s --numpept %s " % (fgfilename.split('/')[-1],bgfilename.split('/')[-1],max(2,numpept))
			

			
	else:
			command1 = "perl c_stat_sig_protein_refpass_depth_memory.pl --foreground %s --background %s" % (fgfilename,bgfilename)
	#print command1
	os.popen(command1)
	
	command2="%sprocess_motifs.pl %slog.%s > %sweboutput%s.p" % (motifPath,motifPath,fgfilename.split('/')[-1]+'.'+bgfilename.split('/')[-1],motifPath,num)
	#print command2
	os.popen(command2)

	outputfile = motifPath+'weboutput%s.p'%num
	try:
		f=open(outputfile,'rb')
		output = f.read()
		f.close()
	except:
		output = ''
	output = output.split('\n')
	data = []
	# sort output on p-value
	for item in output:
		if len(item)>1:
			data.append(  (float(item.split('|')[4].strip()),item))
	data.sort()
## 	########################################

	
	try:
		tablefile = open('%ssubset/motifs/tables/table%s.txt'%(path,num),'w')
	except IOError:
		print "<BR>error openingtablefile<BR>"
		tablefile = open('/tmp/ptmscout/table.txt%s'%random.randint(0,100000000),'w')
	tablefile.write('Motif\tFG\tBG\tP-value\n')

## 	############################################
## 	## display results

	## display query
	try:
		querystr = open("%squery%s"%(outPath,num),'r').read()
		out +=  querystr
	except:
		out += "Error retrieving query string.<br><BR>"
	if len(data)>0:
		out+= '<a href="%ssubset/motifs/foreground.txt?num=%s">Foreground</a>'%(urlpath,num)+'<br>'
		out+= '<a href="%ssubset/motifs/background.txt?num=%s">Background</a>'%(urlpath,num)+'<br>'

		out += '<center><a href="%ssubset/motifs/table.txt?num=%s">Export results to file</a></center>'%(urlpath,num)
		out+= '<center><table   width = 45%% >'
		if type=="peptide":
			out+= '<tr bgcolor="cccccc"><th>Motif</th><th>In Foreground</th><th>Peptides In Background</th><th>P-value</th></tr>'
		else:
			out+= '<tr bgcolor="cccccc"><th>Motif</th><th>In Foreground</th><th>Proteins In Background</th><th>P-value</th></tr>'
			
		for item in data:
			formattedP = "%.2e" % item[0]
			base, power = formattedP.split('e')
			pval = base+' x ' + '10<sup>%s</sup>' % str(int(power))
			realmotif = item[1].split('|')[1]
			motif = realmotif
			inFG = item[1].split('|')[2].split('/')[0]
			inBG = item[1].split('|')[3].split('/')[0]
			out+= '<tr align="center">'
			out+= '<td><tt>%s</tt></td><td>%s</td><td>%s</td><td>%s</td>' %(realmotif.strip(),inFG+'/'+str(numfg),str(inBG)+'/'+str(numbg),pval)
			out+= '</tr>'
			tablefile.write(realmotif.strip()+'\t'+str(inFG)+'/'+str(numfg)+'\t'+str(inBG)+'/'+str(numbg)+'\t'+str(item[0])+'\n')
		out+= '</table></center>'
	else: print "<h3>No results yet. Results may be ready later.</h3>"
	tablefile.close()
	out+= "<BR>"*5		
	fgf = fgfilename.split('/')[-1]
	bgf = bgfilename.split('/')[-1]
	intermediate = 'weboutput%s.p'%num
	intermediate2 = 'log.%s.%s'%(fgf,bgf)
	command4 = "rm -f %s%s" %(motifPath,fgf)
	command5 = "rm -f %s%s" % (motifPath,bgf)
	command6 = "rm -f %s%s" %(motifPath,outputfile)
	command7 = "rm -f %s%s" %(motifPath,fgf+'.'+bgf)
	command8 = "rm -f %s%s" %(motifPath,intermediate)
	command9 = "rm -f %s%s" %(motifPath,intermediate2)
	os.popen(command4)
	os.popen(command5)
	os.popen(command6)
	os.popen(command8)
	os.popen(command9)


	out+= "</body></html>"
	out =  "<html><head><title>PTMScout - Motif Tool</title></head><body>"+out
	try:
		f = open(motifPath+'htmls/'+'outHTML'+str(num),'w')
		f.write(out)
		f.close()
	except IOError: print "outhtml part not working yet"
	########################################
def fifteenmerize(peptide):
	if len(peptide)==15: return peptide
	if 'y' in peptide:
		residue = 'y'
	elif 't' in peptide:
		residue = 't'
	elif 's' in peptide:
		residue = 's'
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
