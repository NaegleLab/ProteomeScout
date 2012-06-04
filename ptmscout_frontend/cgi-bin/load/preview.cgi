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
from template import *
from prevFcns import *
def main():
	form = cgi.FieldStorage()
	filename = form.getvalue("file","")
	try:
		text = open(filename,'r').read()
	except IOError:
		text = ''
		##DEBUG
		filename = "/data/ptmscout/datasets/white_7TimePt_MRM.txt"
		text = open("/data/ptmscout/datasets/white_7TimePt_MRM.txt",'r').read()
	text = text.strip()
	text = text.replace('\r','\n')
	printHeader("PTMScout - Data load preview")
	if text!='':
		try:
			dcols,rcol=getColumns(text)
		except TypeError:
			dcols = []
			rcol = None
		types,labels =  getTypes(dcols,rcol,text)
		testPep =  previewPeptide(text,pepColumn(text))
		pepLines = getTestLines(testPep,text)
		dataDict = getDataToGraph(pepLines,types,labels,rcol,dcols)
		runs,types = getType(text)
		print "<BR><BR><BR><h1>Data graph preview</h1><BR>"
		print helper("Load_Dataset#Preview_Data")+"<BR><BR>"
		print "<table width='500' border='1'>"
		print "<tr valign='top'><td>Runs: </td><td>"
		for item in runs: print "%s<br>"%item
		print "</td></tr>"
		print "<tr valign='top'><td>Data types: </td><td>"
		for item in dcols:
			if 'stddev' not in item[0]:
				print "%s<br>"%item[0]
		print "</td></tr>"
		print "<tr><td>Standard deviations detected</td><td>"
		if "stddev" in text:
			if not (len([item for item in dataDict[dataDict.keys()[0]]['yerrors'] if item == 9])==len(dataDict[dataDict.keys()[0]]['yerrors'])):
				print "YES"
			else: print "NO"
		print "</td></tr>"
		print "</table>"
		if len(dcols)>0:
			print "<img height='600' alt='graph' src = '%simageScripts/preview2.py?file=%s'>"%(urlpath,filename)
		else:
			print "<BR><BR>No data found<BR><BR>"
	else:
		print "<BR><BR>Error loading file."
	print "<BR><BR>Not what you expected to see? See <a href='%stitle=Load_Dataset#Common_Problems' target='_blank'>common problems</a><BR><BR>"%helpPath
	printFooter()

	

main()

