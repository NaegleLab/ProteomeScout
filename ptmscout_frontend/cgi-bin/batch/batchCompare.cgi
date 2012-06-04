#!/usr/local/bin/python
from pathname import *
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import sys

sys.path.append(path+'includes')
from template import *
from tableFunctions import *
def main():
	formdata =''
 	db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    	c=db.cursor()
	form = cgi.FieldStorage()
	expid=form.getvalue("expid",None)
    	printHeader("Batch Comparison Tool", expid)
	

	expTests = [item.replace('include','') for item in form.keys() if not 'id' in item and item != "None"]
	k= len(expTests)
	printExpHeader(expid,c)
	ambiguousPeps = []
	print """<a href="batch.py?exp_id=%s">Back to Comparison Page</a>"""%(expid)
	print """<br><br><A href="%sexperiment.cgi?expid=%s">Back to Experiment Page</a><br><Br>"""%(urlpath,expid)
	
	if k == 0: print "<h2>No experiment selected for comparison</h2>"
	elif k == 1:
		number = 1
		query = """select name from experiment where id = %s"""%expid
		c.execute(query)
		x=c.fetchall()
	 	if len(x)>0: nameExp = x[0][0]
		else: nameExp = ''
		exp=expTests[0]
		exp1=expid
		query = """select name from experiment where id = %s"""%exp
		c.execute(query)
		x=c.fetchall()
		if len(x) >0:
			name = x[0][0]
			exp = (exp,name)
		else: 
			exp = (exp,'')
			name = ''
        	print '<h3>Site Comparison with: '+exp[1]+'</h3><br><br>'
        	exp1Peps = Set(getPeps(exp1,c))
        	exp2Peps = Set(getPeps(exp[0],c))
        	amb = getAmbiguousPeps(exp1Peps,c,exp1)
        	ambiguousPeps=Set([item[0] for item in getAmbiguousPeps(exp1Peps,c,expid)])        	## now add ambiguous peps we know are still for sure in the other db, i.e., if all the possible
        	## resolutions of the ambiguity are also in the other, add to the intersection
        	sd = exp1Peps.symmetric_difference(exp2Peps)
        	
        	for item in amb:
			possible = getPossiblePeps(item[1],c)
			if Set(possible).intersection(Set(exp2Peps))!=None:
            			if item[0] in exp1Peps:
					## if ANY possible ambiguity resolution options in other dataset, remove it from novel sites
					exp1Peps.remove(item[0])
					ambiguousPeps.add(item[0])
            		 
                
                
        	intersect = exp2Peps.intersection(exp1Peps).difference(ambiguousPeps) ## contains all peptides definitely in both
        	exp2only = exp2Peps.difference(intersect) ## contains all peptides only in the big dataset
        	exp1only = exp1Peps.difference(intersect) ## contains all peptides only in the small dataset
        	if len(name)> 30:
            		name = name[0:30]
        	if len(nameExp)>30:
            		nameExp = nameExp[0:30]
		commonAll = exp1Peps
        	if len(sd)!=0:
            		vennNumbers =  [len(exp1only),len(exp2only),len(intersect),len(ambiguousPeps)]
            		print '<img src="%svenn.py?exp1=%s&exp2=%s&intersect=%s&ambiguous=%s&title=%s&number=%s&title2=%s">'% (graphDir,str(vennNumbers[0]),str(vennNumbers[1]),str(vennNumbers[2]),str(vennNumbers[3]),name,str(number),nameExp)
			
			counter = 1
			print "<table width = '50%'>"
			
			for e in expTests:
				query = """select name from experiment where id = %s"""%e
				c.execute(query)
				x=c.fetchall()[0][0]
				exp2Peps = Set(getPeps(e,c))
				common = exp1Peps.intersection(exp2Peps)
				commonAll = commonAll.intersection(exp2Peps)
				print "<tr><td>%s: <a target='_blank' href='common.cgi?exp1=%s&exp2=%s'>%s sites in common</a><br></td></tr>" %(x, expid,e,len(common))
				counter +=1
			### add link for sites common across all experiments
			## make string out of expTests
			expString = ''
			for e in expTests: expString += str(e)+':'
			expString += str(expid)
			#print "<tr><td>All experiments: <a target='_blank' href='common.cgi?expString=%s'>%s sites in common</a><br></td></tr>"%(expString,len(commonAll))
			print "</table>"
			print "<center><h3>Novel sites " + "("+str(len(exp1only))+")"+"</h3></center><br>"
            		makeTable(exp1only,c,exp1,formdata,'ALL',form)
        	else: print "Experiments Identical"
		# regular batch like before
	elif k == 2:
		query = """select name from experiment where id = %s"""%expid
		c.execute(query)
		name1 = c.fetchall()[0][0].replace(' ','%20')
		query = """select name from experiment where id = %s"""%expTests[0]
		c.execute(query)
		name2 = c.fetchall()[0][0].replace(' ','%20')
		query = """select name from experiment where id = %s"""%expTests[1]
		c.execute(query)
		name3 = c.fetchall()[0][0].replace(' ','%20')
		
		exp1Peps=Set(getPeps(expid,c))
		exp2Peps=Set(getPeps(expTests[0],c))
		exp3Peps=Set(getPeps(expTests[1],c))
		intersectionAll = exp1Peps.intersection((exp2Peps).intersection(exp3Peps))
		inAnyIntersection = intersectionAll
		int12 = (exp1Peps.intersection(exp2Peps)).difference(intersectionAll)
		inAnyIntersection.add(int12)
		int13 = (exp1Peps.intersection(exp3Peps)).difference(intersectionAll)
		inAnyIntersection.add(int13)
		int23 = (exp2Peps.intersection(exp3Peps)).difference(intersectionAll)
		inAnyIntersection.add(int23)
		only1 = exp1Peps.difference(inAnyIntersection)
		only2 = exp2Peps.difference(inAnyIntersection)
		only3 = exp3Peps.difference(inAnyIntersection)
		
		print """<img src="%svenn3.py?name1=%s&name2=%s&name3=%s&int12=%s&int23=%s&int13=%s&intall=%s&only1=%s&only2=%s&only3=%s">"""%(graphDir,name1,name2,name3,len(int12),len(int23),len(int13),len(intersectionAll),len(only1),len(only2),len(only3))
		counter = 1
		print "<table width = '50%'>"
		exp1Peps = only1
		ambiguousPeps=Set([item[0] for item in getAmbiguousPeps(exp1Peps,c,expid)])
		commonAll = Set(getPeps(expid,c))
		for e in expTests:
			common = Set(getPeps(expid,c)).intersection(exp2Peps)
			query = """select name from experiment where id = %s"""%e
			c.execute(query)
			x=c.fetchall()[0][0]
			exp2Peps = Set(getPeps(e,c))
			ambiguousPeps=ambiguousPeps.difference(exp2Peps)
			exp1Peps = exp1Peps.difference(exp2Peps)
			amb = getAmbiguousPeps(exp1Peps,c,expid)
        		for item in amb:
            			
            			possible = getPossiblePeps(item[1],c)
				if Set(possible).intersection(Set(exp2Peps))!=None:
            		
					
					if item[0] in exp1Peps:
						## if ANY possible ambiguity resolution options in other dataset, remove it from novel sites
						exp1Peps.remove(item[0])
						ambiguousPeps.add(item[0])
		## amgPeps=Set()
			
			commonAll = commonAll.intersection(exp2Peps)
			print "<tr><td>%s: <a target='_blank' href='common.cgi?exp1=%s&exp2=%s'>%s sites in common</a><br></td></tr>" %(x, expid,e,len(common))
			counter +=1
		expString = ''
		for e in expTests: expString += str(e)+':'
		expString += str(expid)
		print "<tr><td>All experiments: <a target='_blank' href='common.cgi?expString=%s'>%s sites in common</a><br></td></tr>"%(expString,len(commonAll))
		print "</table>"	
			
		
		msids = getMSids_batch(exp1Peps,c,expid)
		newms = []
		for item in msids:
			acc = getTableValues(item, c)[1]
			newms.append(item)
		print '<table width="60%%"><tr><td>'
		
		print '</center><BR>'
		if len(newms)>0:
			print "<center><h3>Novel sites " + "("+str(len(exp1Peps))+")"+"</h3></center><br>"
			makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","pep_aligned"],"",msids=alphabetizeMS(newms,c),width='90',newtab=True)
			
		print '<BR><BR><BR>'
		
		print '</td></tr></table>'
		# figure out how to do venn diagram with 3 things
	else:
		expList=expTests
		exp1Peps = Set(getPeps(expid,c))
		exp1Copy = copy.copy(exp1Peps)
		ambiguousPeps=Set()
		counter = 1
		print "<br><br><table width='60%'>"
		ambiguousPeps=Set([item[0] for item in getAmbiguousPeps(exp1Peps,c,expid)])
		commonAll = exp1Peps
		for ex in expList:
			### print: experiment, number of sites in common, with link to those sites
			exp2Peps = Set(getPeps(ex,c))
			ambiguousPeps=ambiguousPeps.difference(exp2Peps)
		
			query = """select name from experiment where id = %s"""%ex
			c.execute(query)
			x=c.fetchall()[0][0]
			common = exp1Copy.intersection(exp2Peps)
			commonAll = commonAll.intersection(exp2Peps)
			print "<tr><td>Experiment: %s: <a target='_blank' href='common.cgi?exp1=%s&exp2=%s'>%s sites in common</a><br></td></tr>" %(x, expid,ex,len(common))
			
			
			
			amb = getAmbiguousPeps(exp1Peps,c,expid)
			
			exp1Peps = exp1Peps.difference(exp2Peps)
			
	        	
        		
        		for item in amb:
				possible = getPossiblePeps(item[1],c)

				if Set(possible).intersection(Set(exp2Peps))!=None:
            		
					
					if item[0] in exp1Peps:
						## if ANY possible ambiguity resolution options in other dataset, remove it from novel sites
						exp1Peps.remove(item[0])
						ambiguousPeps.add(item[0])
		## amgPeps=Set()
		expString = ''
                for e in expTests:
			expString += str(e)+':'
		expString += str(expid)
		print "<tr><td>All experiments: <a target='_blank' href='common.cgi?expString=%s'>%s sites in common</a><br></td></tr>"%(expString,len(commonAll))
		print "</table><br><br>"
		
        	print "<center><h3>Novel Sites (%s)</h3></center>"%str(len(exp1Peps))
		#print """<a href = '#top'>Top</a>"""
		msids = getMSids_batch(exp1Peps,c,expid)
		newms = []
		for item in msids:
			acc = getTableValues(item, c)[1]
			newms.append(item)
		print '<table width="60%%"><tr><td>'
		
		print '</center><BR>'
		if len(newms)>0:
			makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","pep_aligned"],"",msids=alphabetizeMS(newms,c),width='90',newtab=True)
			
		print '<BR><BR><BR>'
		
		print '</td></tr></table>'
		# do novel site thing with only the experiments in expList
	if len(ambiguousPeps)>0:
			print "<h3>Ambiguous Sites (%s)</h3>"%len(ambiguousPeps)
			print "The following sites have ambiguous protein assignments, some of which may indicate a novel site. <a href='%sambiguity/ambiguity.cgi?expid=%s'><br>Change assignments of ambiguous peptides.</a>"%(urlpath,expid)
			#print "<a href='#top'>Top</a>"
			msids = getMSids_batch(ambiguousPeps,c,expid)
			newms = []
			for item in msids:
				acc = getTableValues(item, c)[1]
		                #print acc.upper()[0], letter
				newms.append(item)
		
			print '<table width="60%%"><tr><td>'
			print '</center><BR>'
			if len(newms)>0:
				makeWholeExpTable(expid,c,"no search","Low",form,["MSid","protein","sequence","site","pep_aligned"],"",msids=alphabetizeMS(newms,c),width='90',icons=False,syns=False,newtab=True)
			print '<BR><BR><BR>'
		
			print '</td></tr></table>'
	printFooter()
 

main()
