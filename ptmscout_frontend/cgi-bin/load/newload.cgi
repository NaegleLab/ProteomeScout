#!/usr/local/bin/python
from pathname import *
from smtplib import SMTP

import sys
sys.path.append(path)
sys.path.append(path+'includes')
from template import *
print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">"""
def main():
	form = cgi.FieldStorage()
	printHeader("PTMScout - Load Data")
	print "<BR><BR>"
	db=MySQLdb.connect(user="%s"%user,passwd="%s"%mysql,db="%s"%database)
	c=db.cursor()
	printFormJS()
	printFillInJS(form,c)
	printKeepFormJS(form)
	
	printForm(c,form)
	print """<script type = 'text/javascript'>fillIN();</script>"""
	printFooter()

def printKeepFormJS(form):
	checkDict = {"true": "true","false":"false"}
	print "<script type= 'text/javascript'>"
	print """
// set the radio button with the given value as being checked
// do nothing if there are no radio buttons
// if the given value does not exist, all the radio buttons
// are reset to unchecked
function setCheckedValue(radioObj, newValue) {
	if(!radioObj)
		return;
	var radioLength = radioObj.length;
	if(radioLength == undefined) {
		radioObj.checked = (radioObj.value == newValue.toString());
		return;
	}
	for(var i = 0; i < radioLength; i++) {
		radioObj[i].checked = false;
		if(radioObj[i].value == newValue.toString()) {
			radioObj[i].checked = true;
		}
	}
}
	function fillIN(){
	    var form = document.forms[0];
	    setCheckedValue(form.loadType,'%s');
	    clearForm(form.loadType,true);
	    form.youremail.value = '%s';
	    form.change.value = '%s';
	    form.link.value = '%s';
	    form.extendChange.value = '%s';
	    form.dataFile.value = '%s';
	    form.description.value = '%s';
	    form.email.value = '%s';
	    form.pmid.value = '%s';
	    form.url.value = '%s';
	    form.published.value = '%s';
	    form.ambiguity.value = '%s';
	    form.modY.checked = %s;
	    form.modT.checked = %s;
	    form.modK.checked = %s;
	    form.modS.checked = %s;
	    form.notes.value = '%s';
	    form.agree.checked = %s;
	    form.loadpm.checked = %s;
	}
	"""%(form.getvalue("loadType","new"),form.getvalue("youremail",""),form.getvalue("change",""),form.getvalue("link",""),form.getvalue("extendChange",""),form.getvalue("dataFile",""),form.getvalue("description",""),form.getvalue("email",""),form.getvalue("pmid",""),form.getvalue("url",""),form.getvalue("published","NO"),form.getvalue("ambiguity","NO"),checkDict.get(form.getvalue("modY","off"),"false"),checkDict.get(form.getvalue("modT","off"),"false"),checkDict.get(form.getvalue("modK","off"),"false"),checkDict.get(form.getvalue("modS",""),"false"),form.getvalue("notes",""),checkDict.get(form.getvalue("agree",""),"false"),form.getvalue("loadpm","false"))
	print "</script>"
def printForm(c,form):
	if form.getvalue("pmidFetch","false") == "true":
		## get author,pub_date,pages,name,journal
		try:
			journal,date,volume,pages,author,name = getPublicationInfo(form,c)
		except ValueError:
			journal,date,volume,pages,author,name= "","","","","",""
			print "<font color='red'>Error retrieving data from pub med</font><BR><BR>"
		except IndexError:
			journal,date,volume,pages,author,name= "","","","",""
			print "<font color='red'>Error retrieving data from pub med</font><BR><BR>"
			
	else:
		journal,date,volume,pages= "","","",""
		author = form.getvalue("authors","")
		name = form.getvalue("expname","")
	
	"""printForm($database_handle,$form) prints the form for loading a dataset, and if there is form data already in the url, it fills it in"""
	print """<table width='60%%'><tr><td>Need help filling out this form? %s</td></tr></table><br><br>"""%helper("Load_Dataset#Dataset_Load_Form")
	print """
	<form action='load.cgi' enctype="multipart/form-data" method="POST" onSubmit="return validateFormOnSubmit(this);">
	<table>
	<tr><td>User Information:</td><td></td></tr>
	<tr><td>Name:</td><td><input type=text name='name' value='%s'></td></tr>
	<tr><td>Email:<font color='red'>*</font></td><td><input type=text name='youremail' value='%s'></td></tr>
	<tr><td><BR><BR>Dataset Information: <font color='red'>*</font></td><td></td></tr>
	<tr><td><br>Type of dataset load %s<BR></td><td></td></tr>
	<tr><td colspan = '2'><br>
	<input onclick="clearForm(this,false);" type=radio value='new' name='loadType' checked>New Dataset<br><BR>
	<input onclick="clearForm(this,false);" type=radio value='reload' name='loadType'>Reload Dataset:<br>%s<br><BR>
	<input onclick="clearForm(this,false);" type=radio value='exportLink' name='loadType'>Extension/change of current experiment:<BR>(If choosing this option you must enter a description of the change in the box below.)<br>%s<BR><BR></td></tr>
	<tr><td>Description of change to extended dataset:</td><td><input type=text name='extendChange' value='%s'> %s</td></tr>
	<tr><td>Data input file:<font color='red'>*</font></td><td><input type=file name='dataFile'> %s <br>Note: files MUST be saved as tab-delimited<BR><BR></td></tr>

		
	<tr><td>Description:<font color='red'>*</font></td><td><textarea rows = '5' cols='50' name='description' >%s</textarea></td></tr>
	<tr><td>Author Contact (email)<font color='red'>*</font></td><td><input type=text size=50 name='email' value='%s'></td></tr>
	<tr><td>PubMed ID:</td><td><input type='text' name='pmid' value='%s'>&nbsp;&nbsp;<input type=checkbox name=loadpm onchange="pmidFetch()">Load citation information automatically.</td></tr>
	<tr><td>URL:</td><td><input type=text name='url' value='%s'></td></tr>
	<tr><td>Published:<font color='red'>*</font></td><td><select name='published'><option value='NO'>NO</option><option value='YES'>YES</option></select></td></tr>
	<tr><td>Experiment Name<font color='red'>*</font></td><td><input type='text' name='expname' value='%s'></td></tr>
	<tr><td>Authors<font color='red'>*</font></td><td><input type='text' name='authors' value='%s'></td></tr>
	<tr><td>Journal<sup>+</sup></td><td><input type=text name = 'journal' value = '%s'></td></tr>
	<tr><td>Publication date<sup>+</sup></td><td><input type='text' name = 'date' value = '%s'></td></tr>
	<tr><td>Volume<sup>+</sup></td><td><input type='text' name ='volume' value = '%s'></td></tr>
	<tr><td>Pages<sup>+</sup></td><td><input type='text' name ='pages' value = '%s'></td></tr>
	<tr><td>+ = required if experiment is published and not loading citation automatically</td><td></td></tr>
	
	<tr><td><br>Is this a mass spec experiment <br>with possibly ambiguous <br>accession assignments?<font color='red'>*</font> %s<br><br></td><td><select name='ambiguity'><option value='NO'>NO</option><option value='YES'>YES</option></select></td></tr>
	<tr><td>Primary residue modification<BR><BR></td><td><input type=checkbox name='modY'>Y
    <input type=checkbox name='modT'>T
    <input type=checkbox name='modS'>S
    <input type=checkbox name='modK'>K<BR><BR></td></tr>
	
	
	<tr  valign='top'><td>Notes: </td><td><textarea rows='5' cols='50'  name='notes'>%s</textarea></td></tr>
	<tr><td><input type=checkbox name = 'agree'> I have read and agree to the <a target='blank' href='%sterms.cgi'>terms of use</a> <font color='red'>*</font></td><td></td></tr>
	<tr><td></td><td><input type=submit value="Submit dataset"></td></tr>
	</table><br><br>
	</form>
	"""% (form.getvalue("name",""),form.getvalue("youremail",""),helper("Load_Dataset#Load_Type"),printExpDropDown("change",c),printExpDropDown("link",c),form.getvalue("extendChange",""),helper("Load_Dataset#Load_Type"),helper("Load_Dataset#Data_Input_File"),form.getvalue("description",""),form.getvalue("email",""),form.getvalue("pmid",""),form.getvalue("url",""),name,author,journal,date,volume,pages,helper("Load_Dataset#Ambiguity"),form.getvalue("notes",""),urlpath)
	



def printFormJS():
	"""printFormJS() prints the validation javascript for the load dataset form"""
	print """<script type="text/javascript">
	// return the value of the radio button that is checked
// return an empty string if none are checked, or
// there are no radio buttons
function getCheckedValue(radioObj) {
	if(!radioObj)
		return "";
	var radioLength = radioObj.length;
	if(radioLength == undefined)
		if(radioObj.checked)
			return radioObj.value;
		else
			return "";
	for(var i = 0; i < radioLength; i++) {
		if(radioObj[i].checked) {
			return radioObj[i].value;
		}
	}
	return "";
}
	function pmidFetch(){
	    var form = document.forms[0];
	    if(form.loadpm.checked==false){return "";}
	    if(form.pmid.value == ""){
	        alert("No pub med ID was entered");
		form.loadpm.checked = false;
	        return "";
		}
	    urlString = "newload.cgi?pmidFetch=true&";
	    for(var i = 0; i < form.elements.length; i++){
	    if (form.elements[i].name != "loadType"){
	       if (form.elements[i].name.indexOf("mod")!=-1 || form.elements[i].name.indexOf("agree")!=-1 || form.elements[i].name.indexOf("loadpm")!=-1){
	          urlString += "&"+form.elements[i].name + "=" + form.elements[i].checked;
	       }
	       else if(form.elements[i].name.indexOf("published")!=-1){
	          urlString += "&"+form.elements[i].name + "="+"YES";
	       }
	       else{
	       urlString += "&"+form.elements[i].name + "=" + form.elements[i].value;
	       }
	       }
	    }
	    urlString += "&loadType=" + getCheckedValue(form.loadType);
	    location.href = urlString;
	}
	function validateFormOnSubmit(theForm) {
	var reason = "";
	reason += validateEmpty(theForm.email.value);
	reason += validateEmpty(theForm.youremail.value);
	reason += validateEmpty(theForm.dataFile.value);
	reason += validateEmpty(theForm.expname.value);
	reason += validateEmpty(theForm.authors.value);
	reason += validateEmpty(theForm.description.value);
	reason+=echeck(theForm.email.value);
	reason+=echeck(theForm.youremail.value);
	reason += checkFile(theForm.dataFile.value);
	reason += validateChecked(theForm.agree.checked);
	reason += validatePublished(theForm);
	if(getCheckedValue(theForm.loadType)=='exportLink'){
	if (theForm.extendChange.value==""){ reason += "\\nMust enter change to the experiment\\n";}
	}
	if (reason != ""){
	
	alert("Some fields need correction:\\n"+reason);
	return false;}
	else
	{return true;}
	}

	function validatePublished(theForm){
	if (theForm.published.value == "YES"){
		if(theForm.pmid.value!=0 && theForm.pmid.value !="" && theForm.loadpm.checked == true){
		return "";
		}
		else if (theForm.pmid.value == "" && theForm.loadpm.checked== true){
		return "Invalid pub med ID";
		}
		else if(theForm.journal.value =="" || theForm.date.value=="" || theForm.pages.value == "" || theForm.volume.value == ""){
		return "Must enter publication data or load from pub med ID";
		}
		else{
		return "";
		}
		
	    }
	    else {
	    return "";
	    }
	}
	function validateChecked(value){
	if (value == true){return "";}
	else {return "\\nYou must agree to the terms of use\\n";}
	}
	function validateEmpty(value){
	if (value == ""){return "Required field left empty\\n";}
	else {return "";}
	
	}
	function checkFile(value)
	{
	 if (value.length>4)
	 {
	 var f = value.split('.').pop();
	 if (f != "txt" && f!="xls" && f!="csv" && f!="ods")
	 {return "\\nInvalid file type:\\n need .txt, ods, .xls, or .cvs\\n";}
	   else {return "";}
	 }
	 else {return "\\nInvalid file type:\\n need .txt, ods, .xls, or .cvs\\n";}
	}

	// return the value of the radio button that is checked
// return an empty string if none are checked, or
// there are no radio buttons
function getCheckedValue(radioObj) {
	if(!radioObj)
		return "";
	var radioLength = radioObj.length;
	if(radioLength == undefined)
		if(radioObj.checked)
			return radioObj.value;
		else
			return "";
	for(var i = 0; i < radioLength; i++) {
		if(radioObj[i].checked) {
			return radioObj[i].value;
		}
	}
	return "";
}

// set the radio button with the given value as being checked
// do nothing if there are no radio buttons
// if the given value does not exist, all the radio buttons
// are reset to unchecked
function setCheckedValue(radioObj, newValue) {
	if(!radioObj)
		return;
	var radioLength = radioObj.length;
	if(radioLength == undefined) {
		radioObj.checked = (radioObj.value == newValue.toString());
		return;
	}
	for(var i = 0; i < radioLength; i++) {
		radioObj[i].checked = false;
		if(radioObj[i].value == newValue.toString()) {
			radioObj[i].checked = true;
		}
	}
}

	

function echeck(str) {

		var at="@"
		var dot="."
		var lat=str.indexOf(at)
		var lstr=str.length
		var ldot=str.indexOf(dot)
		if (str.indexOf(at)==-1){
		   
		   return "Invalid E-mail ID"
		}

		if (str.indexOf(at)==-1 || str.indexOf(at)==0 || str.indexOf(at)==lstr){
		   return "Invalid E-mail ID"
		  
		}

		if (str.indexOf(dot)==-1 || str.indexOf(dot)==0 || str.indexOf(dot)==lstr){
		    return "Invalid E-mail ID"
		    
		}

		 if (str.indexOf(at,(lat+1))!=-1){
		    return "Invalid E-mail ID"
		    
		 }

		 if (str.substring(lat-1,lat)==dot || str.substring(lat+1,lat+2)==dot){
		    return "Invalid E-mail ID"
		    
		 }

		 if (str.indexOf(dot,(lat+2))==-1){
		    return "Invalid E-mail ID"
		   
		 }
		
		 if (str.indexOf(" ")!=-1){
		    return "Invalid E-mail ID"
		    
		 }

 		 return ""
	}


</script>"""

def printExpDropDown(name,c):
	"""printExpDropDown($name,$database_handle) prints an experiment drop down box, with form element name $name"""
	string= """<select onchange="fillForm(this)" name=%s disabled>"""%name
	string+="""<option value='null'>None</option>"""
	query = """select id, name from experiment"""
	c.execute(query)
	x=c.fetchall()
	if len(x)>0:
		for item in x:
			if len(item[1])>100:
				expName=item[1][:100]
			else: expName = item[1]
			string+= "<option value='%s'>%s</option>"%(item[0],expName)
	string+= """</select>"""
	return string
##need dictionary with expid: [expname, authors, author email, description,, pmid, url, published, ambiguous, exportable]

def makeExpDict(c):
	"""makeExpDict($database_handle) makes a dictionary of key=experiment_id, value=another dictionary containing the name, author, contact, description, PMID, URL, published, ambiguity, and export fields for an experiment"""
	query = """select id, name, author, contact, description, PMID, URL, published, ambiguity,export,journal,pages,date,volume from experiment"""
	c.execute(query)
	x=c.fetchall()
	expDict={}
	for item in x:
		id,name,author,contact,description,PMID,URL,published,ambiguity,export,journal,pages,date,volume =item
		d = {}
		d["name"]=name
		d["author"]=author
		d["contact"]=contact
		d["description"]=description
		d["PMID"]=PMID
		d["URL"]=URL
		d["published"]=published
		d["ambiguity"]=ambiguity
		d["export"]=export
		d["journal"]=journal
		d["date"]=date
		d["pages"]=pages
		d["volume"]=str(volume)
		expDict[id]=d
	return expDict

def makeExpJS(expDict):
	"""makeExpJS(expDict) prints the javascript necessary to fill in the form when an experiment is chosen from the drop down box"""
	string = ''
	tfDict = {'0':'NO', '1':'YES'}
	for item in expDict.keys():
		name,author,contact,description,PMID,URL,published,ambiguity,export,journal,date,pages,volume=expDict[item]["name"],expDict[item]["author"],expDict[item]["contact"],expDict[item]["description"],expDict[item]["PMID"],expDict[item]["URL"],expDict[item]["published"],expDict[item]["ambiguity"],expDict[item]["export"],expDict[item]["journal"],expDict[item]["date"],expDict[item]["pages"],expDict[item]["volume"]
		id=item
		string+= """
		if (exp==%s)
		{
		form.expname.value='%s';
		form.authors.value='%s';
		form.email.value='%s';
		form.description.value='%s';
		form.pmid.value='%s';
		form.url.value='%s';
		//form.published.value='';
		//form.ambiguity.value='';
		//form.export.value='';
		form.published.options[%s].selected=true;
		form.ambiguity.options[%s].selected=true;
		//form.export.options[%s].selected=true;
		form.journal.value = '%s';
		form.pages.value='%s';
		form.date.value='%s';
		form.volume.value='%s';
		}
		"""%(id,name,author.replace("'",""),contact,description,PMID,URL,published,ambiguity,export,journal,pages,date.strip(),volume)
	return string
def printFillInJS(form,c):
	"""printFillInJS($database_handle) prints the javascript used to fill in the form when an experiment is chosen"""
	d=makeExpDict(c)
	dictJS = makeExpJS(d)
	print """<script type='text/javascript'>"""
	print"""
	function fillForm(sel)
	{
	var exp = sel.value;
	var form = document.forms[0];
	%s
	}
	"""%dictJS
	print """
	
	function clearForm(radio,loadingpm)
	{
	form = document.forms[0];
	var check = getCheckedValue(radio);
	if (check == "new")
	{
	form.link.value="None";
	form.link.disabled = true;
	form.change.value="None";
	form.change.disabled=true;
	}
	else if (check == "reload")
	{
	form.link.value="None";
	form.link.disabled = true;
	form.change.disabled=false;
	}
	else
	{
	form.link.disabled = false;
	form.change.value="None";
	form.change.disabled=true;
	}

	if(loadingpm == false){
	form.dataFile.value='';
	
	   form.expname.value='';
	   form.authors.value='';
	   form.journal.value = '';
	   form.volume.value= '';
	   form.pages.value='';
	   form.date.value = '';
	   
	   form.journal.value = '';
	   form.pages.value ='';
	   form.date.value='';
	   form.volume.value='';
        form.email.value='';
	form.description.value='';
	form.pmid.value='';
	form.url.value='';
	form.notes.value='';
	}
	}
	"""
	print """</script>"""

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
		author = lines[0].split("\t")[1]
		name = lines[4].split("\t")[1]
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
	return journal,date,volume,pages,author,name
	
main()
