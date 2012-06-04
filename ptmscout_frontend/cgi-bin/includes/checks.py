import sys
from pathname import *
sys.path.append(path+'imageScripts')
from graph import *

def printJS():
	printAddRemoveFunction()  
	printChangeChecks()

def printCols(cols,name):
	string = ''
        string+= '<table width = 100%%>'
        counter = 2
        string+= '<tr>'
	string+= "<td><input type=checkbox name = '%s' value='all'>all" % name
        for item in cols[1:]:
                string+= '<td>'
                string+= "<input type=checkbox name = '%s' value = '%s'>%s" % (name,item,item)
                string+= '<\/td>'
                if counter %2==0: string+= '<tr>'
                counter +=1
	string+='<\/tr>'
        string+= '<\/table>'
	return string
#1 is mouse, 2 is expression human, 3 is nci60
def printGOTable(d, aspect): # aspect is either 'B','C','F'
	# print table of GO terms - functions table
	string=''
	string+="""<table>"""
	for item in d["go"]:
		# color code according to which GO aspect
		if item[0]==aspect:
			string+= "<tr><td>"
			string+= """<a href='%s""" % url_dict.get('GO')
			string+= """%s' target='_blank'>""" % item[1]
			string+= """<font color= %s>""" % go_dict.get(item[0])
			string+= """ %s</font></a>""" % item[1]
			string+= """ (%s)""" % item[2]
			string+= "</td></tr>"
			string+= '</table>'
	return string

# divs will be of the form aspect initial, protein id, i.e. B9378
def printAddRemoveGORow(d):
	print '<script type="text/javascript">'
	print """
function addRemoveGORow(form,div1,div2,div3)
{
        if (form.downup.value == "show")
	{
	document.getElementById(div1).innerHTML = " """ + printGOTable(d,'F') + """ ";
        document.getElementById(div2).innerHTML = " """ + printGOTable(d,'C') + """ ";
	document.getElementById(div3).innerHTML = " """ + printGOTable(d,'B') + """ ";
	form.downup.value = "hide";
       	}
	else
	{
	document.getElementById(div1).innerHTML = "";
	document.getElementById(div2).innerHTML = "";
	document.getElementById(div3).innerHTML = "";
	form.downup.value = "show";
	}

}
"""
	
def printAddRemoveFunction():
	print '<script type="text/javascript">'
	print """
function addRemoveChecks(form)
{
	if (form.buttons.value== "show column choices")
		{
			if (form.table.value=="null")
				{
				document.getElementById('checks').innerHTML =" """+ str(printCols(EXPRESSION_MOUSE_COLS,"mouse"))+""" ";
				}
			if (form.table.value=="tissue")
				{
				document.getElementById('checks').innerHTML = " """+str(printCols(EXPRESSION_HUMAN_COLS,"human_tissue_cols"))+""" ";
				}
			if (form.table.value=="cell")
				{
				document.getElementById('checks').innerHTML = " """+str(printCols(EXPRESSION_NCI60_COLS,"nci60_cols"))+""" ";
				}
			if (form.table.value=="both")
				{
				document.getElementById('checks').innerHTML = " """+str(printCols(EXPRESSION_NCI60_COLS,"nci60_cols"))+"<br><hr><br>"+str(printCols(EXPRESSION_HUMAN_COLS,"human_tissue_cols"))+""" ";
				}
			form.buttons.value= "hide column choices";
		}
	else
		{
		document.getElementById('checks').innerHTML = "";
		form.buttons.value ="show column choices";
		}
}"""
	print '</script>'

def printChangeChecks():
	print '<script type="text/javascript">'
	print """
function changeChecks(form)
{
	if (form.buttons.value== "hide column choices")
		{
			if (form.table.value=="null")
				{
				document.getElementById('checks').innerHTML =" """+ str(printCols(EXPRESSION_MOUSE_COLS,"mouse"))+""" ";
				}
			if (form.table.value=="tissue")
				{
				document.getElementById('checks').innerHTML = " """+str(printCols(EXPRESSION_HUMAN_COLS,"human_tissue_cols"))+""" ";
				}
			if (form.table.value=="cell")
				{
				document.getElementById('checks').innerHTML = " """+str(printCols(EXPRESSION_NCI60_COLS,"nci60_cols"))+""" ";
				}
			if (form.table.value=="both")
				{
				document.getElementById('checks').innerHTML = " """+str(printCols(EXPRESSION_NCI60_COLS,"nci60_cols"))+"<br><hr><br>"+str(printCols(EXPRESSION_HUMAN_COLS,"human_tissue_cols"))+""" ";
				}
			
		}
	else
		{
		document.getElementById('checks').innerHTML = "";
		
		}
}"""
	print '</script>'

