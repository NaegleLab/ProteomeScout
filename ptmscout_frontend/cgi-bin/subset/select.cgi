#!/usr/local/bin/python
from pathname import *
import sets
import cgi

import cgitb
if displayPythonErrors:
    cgitb.enable()

print "Content-Type: text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import MySQLdb
import sys

sys.path.append(path+'includes')
sys.path.append(path)
sys.path.append(path+'clusters')
from template import *
import pickle
import os



from pathname import *
sys.path.append('clusters')

cluster_file_path=""
user_file_name=""
random_file_key=""

from clusters import *
NUM_DATA_DIVS = 7
MAX_NUM_SUBSETS = 30
def main():
    global cluster_file_path, random_file_key, user_file_name
    
    form = cgi.FieldStorage()
    if form.has_key("exp_id") and form["exp_id"].value!="":
        exp_id=form.getvalue("exp_id")
    else:
        exp_id = '7'
    clusters = form.getvalue("clusters","no")
    printHeader("Evaluate Subsets",exp_id)
    
    db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    c=db.cursor()
    query = """select id from experiment where export = 0"""
    c.execute(query)
    x=c.fetchall()
    exps = [str(item[0]) for item in x]
    if exp_id in exps:
        print "<BR><div class='error'>Illegal experiment id. Subset selection not allowed for this experiment.</div><br><BR><BR>"
        printFooter()
        sys.exit(0)
    print """<BR>%s<BR>"""%(helper("title=Subset_Selection"))
    printExpHeader(exp_id,c)
    
    query = """select id from experiment where export = 0"""
    c.execute(query)
    x=c.fetchall()
    exps = [str(item[0]) for item in x]
        
    ##load cluster file
    ## write it to the clusterPath directory with a random name.  Save the random name
    ## as a variable in this python script. IF we need to refer back to the file later, we
    ## Use the random name (regenerated) instead of the real file path.
    load = form.getvalue('load','no')
    if load == "yes":
        file = form.getvalue('clusterFile','')
        file_item = form['clusterFile']
        if file =="":
            clusters="no"
            print "<div class='error'>No file uploaded.</div>"
        else:
            if checkFile(file) == False:
                print "<div class='error'>Not a valid cluster file.</div>"
            else:
                # Move the uploaded file into the clusters directory
                random_file_key = random.randint(0,10000000)

                # Better here would be to make sure the file doesn't exist before overriding it
                cluster_file_path = hackedFileName(exp_id, random_file_key)

                user_file_name=file_item.filename
                open(cluster_file_path,'wb').write(file)

                if ClusterSet(cluster_file_path).parseError != None:
                    clusters = "no"
                    print "<div class='error'>Error parsing cluster file. See the help page for advice on importing clusters.</div>"
                else:
                    print "<div class='feedback'>Cluster file %s loaded successfully for experiment %s</div>" % (str(user_file_name), str(exp_id))
                    clusters = "yes"
    
    ### end load cluster file


    c.execute("""select * from data join MS on data.MS_id=MS.id where experiment_id=%s""",(exp_id,))
    x=c.fetchall()
    data = True
    if not(len(x)>0): data = False
    printJavaScript(exp_id,c,data,clusters=clusters)
    ## keep hidden counter so we know which div of subsets is activated
    print '<input type = hidden name = counter value = 0>'
    
    
    ### load clusters
    print """
    <h3>Clustering %s</h3>""" %(helper("Subset_Selection#Export_File_For_Clustering"))
    print '<table id="cluster_table" cellspacing="10px" cellpadding="4px">'
    print '<tr>'
    if str(exp_id) not in exps:
        runs = getType(exp_id)[0]
        print """
        <td valign="top">
        <h4>Experiment Export</h4><form method = "POST" enctype="multipart/form-data" action ="clusters/exportFile.txt">
        <input type = hidden name = expid value = %s>
        <input type = submit value = "Export Experiment for Clustering"><BR>"""% (str(exp_id))
        if len(runs)>1:
            print """
            <input type=checkbox name = 'avg'>Average multiple runs
            <br></form>
            *Selecting "average" gives a data column with the average value of all runs. By default, for datasets with mutliple runs, data columns for each run will be given in random order. You may need to change the order of the columns.
            </td>
            """
        else:
            print """
            <br></form></td>"""
                
    print """
    <td valign="top">
    <h4>Cluster Upload</h4>
    <form method = "POST" enctype="multipart/form-data" action = "select.cgi?exp_id=%s">
    """ % (str(exp_id))        
    print """
    <p>
    Cluster File:
    <input type="file" name="clusterFile"> %s
    </p>
    <p>
    """%(helper("Subset_Selection#Export_File_For_Clustering"))
    print """
    <input type = "submit" value = "Load Clusters" /> (Resets current search)


    </p>
    <input type="hidden" name = "exp_id" value = "%s">
    <input type = "hidden" name = "load" value = "yes"></form>""" %(str(exp_id))

    print "</td>"
    
    # okay, here is where we drop them into mcam
    # we create a new HTML FORM element with a bunch of input and a submit button
    # This gets posted back to the server (a different script) where we invoke the
    # perl script to do the processing.

    if clusters=="yes":
        print "<td valign='top'>"
        print "  <h4>Multiple Clustering Analysis Methodology</h4>"
        print """  <form method="POST" action="mcam.cgi">"""
        print """  <input type="hidden" name="file_key" value="%d">"""%(random_file_key)
        print """  <input type="hidden" name="exp_id" value="%s">"""%(str(exp_id))
        print """  <input type="hidden" name="user_file_name" value="%s">"""%(str(user_file_name))
        print "  <input type='Submit' value='Launch MCAM' />"
        print "  </form>"
        print "</td>"
        

    print '</tr></table>'
    ### end load clusters

    ## begin subset selection form ##
    print '<div id="content3">'
    printForm(exp_id,c,data,form,clusters=clusters)
    print '<br><br></div>'
    ## end subset selection form ##

    ## keep hidden list (separated by colons) of msids in currently selected subset
   
    for i in range(1,MAX_NUM_SUBSETS):
        print '<input type = hidden name = selectedids%s value = "">'%i

    ## allow to export tabs as cluster set
    print "<h3>Explore Subsets %s</h3>"%(helper("Subset_Selection#Explore_Subset"))
    print "<input  onClick='Features();'type='button' value='Evaluate features across subsets'>"
    print "<input  onClick='exportClusters();'type='button' value='Export existing tabs as cluster set'><br>"

    ## print tabs (later will be pretty images)
    print '<div id = tabs align = left>'
    print '<table cellspacing = 0><tr>'
    for i in range(1,MAX_NUM_SUBSETS):
        if i %13==0: print "</tr><table><table cellspacing = 0><tr>"
        print '<td>'
        printImageMap(i)
        print '</td>'
        print """<td><div id="font%s" style="position: relative;display:none">
        <p style="position: absolute; margin: 5px; bottom: 25px; top: 0px; left: 10px; padding: 0px; color: white;"><b><a href="javascript:changeDivs(%s)">
        <font size='1' color=white>Subset %s</font></a></b></p><img alt="tab" usemap="\#tabmap%s" width='100' hspace=0 border = none src = "%stab_active.jpg" style="display:none" name = b%s></div></td>""" % (i,i,i,i,imgPath,i)
        
    print '</tr></table>'

    print '</div>'
    subsetFrame(MAX_NUM_SUBSETS)
    db.close()
    ### did we click on a term from summary?
    url = form.getvalue("url",None)
    if url !=None:
        url= url.replace('code1','&')
        url = url.replace('code2','?')
        print """
        <script type=text/javascript>
        var url = '%s';
        document.getElementsByName('sub1')[0].src=url;
        set();
        </script>
        """%url
    ###
    printFooter()

def printForm(exp_id,c,data,form,clusters='no'):
    numMS,numPro,numPep = getExperimentParameters(c,exp_id)
    print "<h3>Subset Selection %s</h3>" %(helper("Subset_Selection#Definitions_in_Subset_Select_Menu"))
    print '<div class="nested_content">'
    print '<span class="help">'
    print 'Use the search form below to perform the following tasks:'
    print '<ul>'
    print '<li>Explore a smaller portion of a large dataset based on a common feature</li>'
    print '<li>Generate hypotheses based on enrichment between dynamics and annotations or annotations of different sources</li>'
    print '<li>View enrichment in clustering imported from outside (e.g. Self organizing maps, kmeans etc.) </li>'
    print '</ul>'
    print '</div>'
    
    print '<div id="content2" class="nested_content subset_search">'
    print '<p class="search_reset">'
    print '<a href="javascript:window.location=document.location">Reset Search</a> |'
    print '<a href="javascript:clearForm();">Clear Form</a>'
    print '</p>'
    print '<p>'
    print """
     

            """
    print """<form  name = form1 target="sub1" method="POST" action ="search.cgi?exp_id=%s" >""" % (exp_id)
    print '<input type="hidden" name = "numMS" value = "%s">' % numMS
    print '<input type="hidden" name = "numPro" value = "%s">' % numPro
    print '<input type="hidden" name = "numPep" value = "%s">' % numPep
    if (random_file_key != ""):
        print '<input type="hidden" name = "file_key" value = "%d">' % int(random_file_key)
        
    print '<input type="hidden" name = "user_file_name" value = "%s">' % str(user_file_name)
    
    print '<div>'
    print """
    
    <fieldset><legend><strong>Filter</strong></legend>

    <p>
    <label for="select">Background:</label>    
    <select name = select >
      <option value = all>Entire Experiment</option>
      <option value = current>Last viewed subset</option>    
    </select>
    </p>
    
    <input type=hidden name=exp_id value = %s>"""%(exp_id)
    print """<input type=hidden name=selectedids value =''>"""

    
    if data:
        print """
        <p>
        <label for="dataType">Criteria Type:</label>
        <select name = dataType onChange = "filterDivs(this.form)">
          <option value = both>All</option>
          <option value = data>Quantitative Measurements</option>
          <option value = metadata>Metadata</option>
        </select>
        </p>
          """
    print"""
     
    </fieldset>
    
    </div>
    
    """
    
    # begin data div
    if(data):
        print '<div id = "data">'
    else:
        print '<div id="data" style="display:none">'
        
    print '<fieldset><legend><strong>Quantitative Measurement Search</strong></legend>'

    dataDivCount = 0
    vals = getValueOptions(exp_id,c,data)
    print makeNewDataDiv('FIRST',dataDivCount,vals)
    dataDivCount = 1
    for i in range(NUM_DATA_DIVS-1):
        print makeNewDataDiv('AND',dataDivCount,vals)
        print makeNewDataDiv('OR',dataDivCount, vals)
        dataDivCount += 1
    print '</fieldset>'
    print '</div>'

    # end data div
    
    print '<div align=left></div>'
    # begin metadata div
    print '<div id = "metadata">'

    print '<fieldset><legend><strong>Metadata Search</strong></legend>'
    print '<p>'
    print '<label for="strin">Scansite Predictions Stringency:</label><BR>'
    print   "<input type = radio name='strin' value='low'>Low Stringency (score &le; 5)<br/>"
    print   "<input type = radio name='strin' value='medium' checked>Medium Stringency (score &le; 1)<br/>"
    print   "<input type = radio name='strin' value='high'>High Stringency (score &le; 0.2) <br/>"
    print "</p>"

    print '<strong>Metadata:</strong>'
    dataDivCount = 0
    print makeNewMetadataDiv('MFIRST',dataDivCount,exp_id,c,clusters = clusters)
    dataDivCount =1
    for i in range(NUM_DATA_DIVS-1):
        print makeNewMetadataDiv('MAND',dataDivCount, exp_id,c)
        print makeNewMetadataDiv('MOR',dataDivCount, exp_id,c)
        dataDivCount += 1
    print '</fieldset>'
    
    print '</div>'
    # end metadata div

    print '<div class="nested_content" style="margin-left:auto; margin-right:auto">'
    print '<input type="submit" value="Search" onclick = "set()">'
    print '</div>'
    print """</form>"""
    print '</div>'
    
# type is either FIRST, AND, or OR (first for the first instance of the div, AND or OR depending on which button they are selected by)
# Each div has: ( I will call value/operator/comparator opType)
# value  operator   value  comparator value (--> operator  value)
# where operator = *, /, +, - and comparator = <, > , = , !=, <=, >=
# each div is numbered as type_number, so the first div is FIRST0
# each form element has name of the form: divid_orderInDiv_opType, so the value
# in the first div has name = FIRST0_0_value

# all values start as drop downs, then can be toggled to become text input fields to input a numerical value
# so ids of values are FIRST0_0_value_div for its container, FIRST0_0_value_dd for the drop down div, FIRST0_0_value_num
# for its text box number field
def getMetaValueOptions(exp_id,c,clusters='no',type="MFIRST"):
    string = ''
    string += '<option value = "None">Null</option>'
    string += '<option value = site_type>Site Type</option>'
    string += '<option value = GO_BP>GO - BP</option>'
    string += '<option value = GO_MF>GO - MF</option>'
    string += '<option value = GO_CC>GO - CC</option>'
    string += '<option value = acc_gene>Protein Name</option>'
    string += '<option value = swissprot>SwissProt ID</option>'
    
    string += '<option value = gi>gi</option>'
    string += '<option value = entrez_protein>entrez</option>'
    string += '<option value = pfam_site>Pfam Site</option>'
    string += '<option value = species>Species</option>'
    string += '<option value = domain>Domain</option>'
    if type == "MFIRST": string += '<option value = pep_aligned>Aligned Sequence</option>'
    string += '<option value = scansite_kinase>Scansite-Kinase</option>'
    string += '<option value = scansite_bind>Scansite-Bind</option>'
    string +='<option value=pelm_kinase>PhosphoElm Kinase</option>'
    if clusters == 'yes' and type == "MFIRST":
        string+='<option value = "clusters">Clusters</option>'
    return string
    
def getValueOptions(exp_id,c,data):
    if not(data): return '<option value = none></option>'
    c.execute("""select run, type, label from data join MS on data.MS_id = MS.id where experiment_id = %s""",(exp_id))
    x=c.fetchall()
    runs = [item[0] for item in x]
    types = [item[1] for item in x]
    labels = [item[2] for item in x]
    #remove duplicates from each list
    runs = reduce(lambda l, x: x not in l and l.append(x) or l, runs, [])
    types = reduce(lambda l, x: x not in l and l.append(x) or l, types, [])
    labels = reduce(lambda l, x: x not in l and l.append(x) or l, labels, [])
    if len(runs) == 1:
        if 'time' in runs[0]:
            druns=[(runs[i],runs[i]) for i in range(len(runs))]
        else: druns=[(runs[i],'') for i in range(len(runs))]
    else: druns = [(runs[i],runs[i]) for i in range(len(runs))]
    if len(types)==1:
        if 'time' in types[0]:
            dtypes=[(item,item) for item in types]
        else:
            dtypes=[(item,'') for item in types]
    else: dtypes = [(item,item) for item in types]
    
        
    # now make every possible combination
    combinations = []
    dcombinations = []
    for run in druns:
        for type in dtypes:
            for label in labels:
                dcom = ''
                
                if run[1] !='':
                    dcom += run[1] + ':'
                if type[1] != '': dcom += type[1]+':'
                dcom += label
                
                combinations.append(type[0]+':'+run[0]+':'+label)
                dcombinations.append(dcom)
    string = ''
    for i in range(len(combinations)):
        string += '<option value = "%s">%s</option>' % (combinations[i], dcombinations[i])
    return string
        
def makeValue(type, order, order2,valOpts):
    string = ''
    string += """<select name = %s onChange = "numBox('%s','%s')">""" % (type+str(order)+'_'+str(order2)+'value',getId(type,order)+str(order2)+'_num',type+str(order)+'_'+str(order2)+'value')
    string += '<option value = "None">Null</option>'
    string += valOpts
    string += '<option value = "num">Numerical</option>'
    
    string += '</select>'
    # end dd for value 0
    return string
def makeOp(type,order,order2):
    string = ''
    string += ''
    string += '<select name = %s>' % (type+str(order)+'_'+str(order2)+'operator')
    string += '<option value = "None">Null</option>'
    string +=  '<option value = "plus">+</option>'
    string +=  '<option value = "minus">-</option>'
    string +=  '<option value = "mult">&times;</option>'
    string +=  '<option value = "divide">&divide;</option>'
    if order2 == 1:
        string+='<option value="gt">></option>'
        string+='<option value="lt"><</option>'
        string+='<option value="geq">&ge;</option>'
        string+='<option value="leq">&le;</option>'
    string += '</select>'
    return string
def makeComp(type,order,order2):
    string = ''
    string += '<select name = %s>' % (type+str(order)+'_'+str(order2)+'comparator')
    string += '<option value = "None">Null</option>'
    string +=  '<option value = "equals">=</option>'
    string +=  '<option value = "notequals">&ne;</option>'
    string +=  '<option value = "gt">></option>'
    string +=  '<option value = "lt"><</option>'
    string +=  '<option value = "geq">&ge;</option>'
    string +=  '<option value = "leq">&le;</option>'
    string += '</select>'
    return string
def getId(type,order):
    return  type + str(order)
def makeNewDataDiv(type,order,vals):
    string = ""
    if type == "FIRST":
        string += '<div class="data_search nested_content" id="%s">' % getId(type,order)
    else:
        string += '<div class="data_search nested_content" id = "%s" style="display:none">' % getId(type,order)
    if type != "FIRST" and order != NUM_DATA_DIVS-1:
        string +=  type
    if type !="FIRST": string += '<input name = %s type = hidden value = "inactive">'% getId(type,order)
    string += '<table><tr>'
    # value 0
    string += '<td>'
    string += makeValue(type,order,0,vals)
    string += '</td>'    
    # operator 1
    string += '<td>'
    string += makeOp(type,order,1)
    string += '</td>'
    # value 2
    string += '<td>'
    string += makeValue(type,order,2,vals)
    string += '</td>'
    # comparator 3
    string += '<td>'
    string += makeComp(type,order,3)
    string += '</td>'
    # value 4
    string += '<td>'
    string += makeValue(type,order,4,vals)
    string += '</td>'
    # operator 5
    string += '<td>'
    string += makeOp(type,order,5)
    string += '</td>'
    # value 6
    string += '<td>'
    string += makeValue(type,order,6,vals)
    string += '</td>'
    string += '</tr>'

    # make possibility to enter numerical value
    string += '<tr><td><div id = "%s" style="display:none"><input type=text name = %s></div></td>' % (getId(type,order)+'0_num',type+str(order)+'_'+str(0)+'value_num')
    string += '<td></td>'
    string += '<td><div id = %s style="display:none"><input type=text name = %s></div></td>' % (getId(type,order)+'2_num',type+str(order)+'_'+str(2)+'value_num')
    string += '<td></td>'
    string += '<td><div id = %s style="display:none"><input type=text name = %s></div></td>' % (getId(type,order)+'4_num',type+str(order)+'_'+str(4)+'value_num')
    string += '<td></td>'
    string += '<td><div id = %s style="display:none"><input type=text name = %s></div></td>' % (getId(type,order)+'6_num',type+str(order)+'_'+str(6)+'value_num')
    string += '</tr>'

    string += '</table>'
    if order != NUM_DATA_DIVS-1:
        if "FIRST" in type:
            string += """<input type=button value = "AND &darr;" onClick="addDataDiv(%s,'AND')">&nbsp;&nbsp; <input type = button value = "OR &darr;" onClick = "addDataDiv(%s,'OR')">&nbsp;&nbsp;"""%(order+1,order+1)
        elif "OR" in type:
            string += """<input type = button value = "OR &darr;" onClick = "addDataDiv(%s,'OR')">&nbsp;&nbsp;"""%(order+1)
        elif "AND" in type:
            string += """<input type=button value = "AND &darr;" onClick="addDataDiv(%s,'AND')">"""%(order+1)
    if type != "FIRST":
        string += '<input type=button value = "CLOSE &uarr;" onClick = "noDataDiv(%s)">' % (order)
    string += '</div>'
    return string


def makeNewMetadataDiv(type,order,exp_id,c,clusters='no'):
    string = ""
    if type == "MFIRST":
        string += '<div class="metadata_search nested_content" id="%s">' % getId(type,order)
    else:
        string += '<div class="metadata_search nested_content" id = "%s" style="display:none">' % getId(type,order)
    if type != "MFIRST" and order != NUM_DATA_DIVS-1:
        string +=  type[1:]
    if type !="MFIRST":
        string += '<input name = %s type = hidden value = "inactive">'% getId(type,order)
   
    string += '<table><tr><td>'
    string += """<select name = %s onChange = "metaDataBox(%s,%s,this.form)">""" % (getId(type,order),getId(type,order)+'_value',getId(type,order))
    string += getMetaValueOptions(exp_id,c,clusters=clusters,type=type)
    string += '</select>&nbsp;&nbsp;'
    string +='</td><td>'
    if type=="MFIRST":
        string += """<select name = %s onChange = "Cluster(this);"><option></option></select>""" % (getId(type,order)+'_value')
        string += """<select name = "cluster%s"  style="display:none;"><option></option></select>"""%(order)
        string += """<select  name = "seq%s" style="display:none;" ><option></option></select>""" % order
    else:
        string += """<select name = %s ><option></option></select>""" % (getId(type,order)+'_value')
    string += '<input type = text name = "pro%s" style="display:none;">' % order
    string += '<select name="stringency" style="display:none"><option value="medium">Medium stringency &le; 1</option><option value="high">High stringency &le; 0.2</option><option value="low">Low stringency &le; 5</option></select>'
    string += '</td></tr></table>'
    string += '<br>'
    if order != NUM_DATA_DIVS-1:
        if "FIRST" in type:
            string += """<input type=button value = "AND &darr;" onClick="addMetaDataDiv(%s,'MAND')">&nbsp;&nbsp; <input type = button value = "OR &darr;" onClick = "addMetaDataDiv(%s,'MOR')">&nbsp;&nbsp;"""%(order+1,order+1)
        elif "AND" in type:
            string += """<input type=button value = "AND &darr;" onClick="addMetaDataDiv(%s,'MAND')">"""%(order+1)
        elif "OR" in type:
            string += """<input type = button value = "OR &darr;" onClick = "addMetaDataDiv(%s,'MOR')">&nbsp;&nbsp;"""%(order+1)
    if type != "MFIRST":
        string += '<input type=button value = "CLOSE &uarr;" onClick = "noMetaDataDiv(%s)"><BR><BR>' % (order)
    string += '<br><br></div>'
    return string






def makeClusterOptionsDict(clusterDict):
    d = {}
    for item in clusterDict.keys():
        if item != '':
            keys = clusterDict.get(item,{}).keys()
            if False not in [not(key.isalpha()) and key.isalnum() for key in keys]:
                keys = [int(key) for key in keys]
            keys.sort()
            d[item] = printClusterOptionBox(keys)
    return d

def makeOptionsDict(exp_id,c,clusters='no'):
    d = dict()
    GO = getGO(exp_id,c)
    d['GO_BP'] = printSelectOptionBox(GO[0],GO=True,addNull=True)
    d['GO_MF'] = printSelectOptionBox(GO[1],GO=True,addNull=True)
    d['GO_CC'] = printSelectOptionBox(GO[2],GO=True,addNull=True)
    d['acc_gene'] = printSelectOptionBox(getAcc(exp_id,c))
    d['swissprot'] = printSelectOptionBox(getSwiss(exp_id,c))
    d['entrez_protein'] = printSelectOptionBox(getEntrez(exp_id,c))
    d['gi'] = printSelectOptionBox(getGi(exp_id,c))
    d['species'] = printSelectOptionBox(getSpecies(exp_id,c))
    d['pfam_site']=printSelectOptionBox(getPfam(exp_id,c))
    d['domain']=printSelectOptionBox(getDomain(exp_id,c))
    d['site_type'] = printSelectOptionBox(['Y','T','S','K'])
    d['scansite_bind']=printSelectOptionBox(getBind(exp_id,c),addNull=True)
    d['scansite']=printSelectOptionBox(getScan(exp_id,c),addNull = True)
    d['scansite_kinase']=printSelectOptionBox(getSKinase(exp_id,c),addNull=True)
    d['pelm']=printSelectOptionBox(getPELM(exp_id,c),addNull=True)
    if clusters=='yes':
        d['clusters']=printSelectOptionBox(getClusters(cluster_file_path))
    else:
        d['clusters']=printSelectOptionBox([])
    return d

def getClusters(filename):
    clusterSets = ClusterSet(filename)
    names = clusterSets.clusterSets.keys()
    if False not in [item.isalnum() and not(item.isalpha()) for item in names]:
        names = [int(item) for item in names]
    names.sort()
    names = ["Cluster set: "+str(item) for item in names]
    return names
   

# each subset div can be one of:
# nonexistent: before a subset is chosen - empy tab
# active: it is the tab in focus - dark tag
# inactive: it has been selected but another tab is being focused on

# on load (put in the onclick of whatever takes you to a new frame) of any frame, set all other existent divs to inactive
# set hidden "currently selected" msids to whatever is in the frame

def subsetFrame(MAX_NUM_SUBSETS):
    for i in range(1,MAX_NUM_SUBSETS+1):
        print '<div id = "subset%s"  style="display:none">' % i
        print '<input type= hidden name = level%s value = "nonexistent">'%i
        print '<iframe  width = "95%%"  frameborder=0 scrolling = "no" src="" name = "sub%i"></iframe>' % (i) 
        print '</div>'


 
def printJavaScript(exp_id,c,data,clusters='no'):
    optionsDict = makeOptionsDict(exp_id,c,clusters=clusters)
    
    if data:
        print '<script type = "text/javascript">'
   
        print """
    function numBox(box,sel)
    {
        if (document.getElementsByName(sel)[0].value == "num")
        {document.getElementById(box).style.display = 'block';} 
        else {document.getElementById(box).style.display='none';}

    }"""
        print '</script>'
    
    # box is name of value box, sel is name of select type box
    print "<script type = 'text/javascript'>"
    print """function exportClusters()
    {
    
    var url = 'clusters/exportClusters.txt?';
    for(i = 1; i<20; i++)
       {
       
       
       
       if (document.getElementsByName("level" + i)[0].value != 'nonexistent')
       {
       var ids = document.getElementsByName('sub'+i)[0].contentWindow.document.msform.msids.value;
       url += 'ids'+i+'='+ids+'&';
       }
       }
       window.open(url);

    }
    """
    print """function Features()
    {
    var url = 'features.py?';
    for(i = 1; i<20; i++)
       {
       
       
       
       if (document.getElementsByName("level" + i)[0].value != 'nonexistent')
       {
       if (document.getElementsByName('sub'+i)[0].contentWindow.document.getElementsByName('featureDict')[0] !=undefined)
       {
       var file = document.getElementsByName('sub'+i)[0].contentWindow.document.getElementsByName('featureDict')[0].value;

       url += 'file'+i+'='+file+'&';
       }
       }
       }
       window.open(url);
       //alert(url);

    }
    """
    print "</script>"
    print '<script type = "text/javascript">'
    print """function Cluster(box){clusterBox('cluster0',box.value);}
    """

    print """function MCAM()
    {
    var url = 'MCAM.py?';
    window.open(url);
    }"""
    
    print '</script>'
    
    if clusters=='yes':
         print '<script type = "text/javascript">'
         clusterSets = ClusterSet(cluster_file_path)
         clusterDict = clusterSets.clusterSets
         
         clusterOptionsDict = makeClusterOptionsDict(clusterDict)
         print """
          function clusterBox(box,key)
          {
          
          if(typeof(box)=='string'){var box = document.getElementsByName(box)[0];}
          if (typeof(key) != 'string'){key = key.value;}
          key = key.replace('Cluster set: ','');
          box.options.length = 0;"""
         print """box.options.add(new Option('All clusters','ALL'));
          """
         
         for item in clusterOptionsDict.keys():
             
             
             print "if(key=='%s'){%s}" % (item.strip(),clusterOptionsDict[item])
             

         print """
           }

           """
         print '</script>'
    
    else: print "<script type='text/javascript'>function clusterBox(box, sel){var i = 0;}</script>"
    print '<script type = "text/javascript">'
    print """
    function metaDataBox(box,sel,form)
    {
        
        
        var string1 = box.name;
        var n = '';
        
        for(i = 0; i < string1.length; i++)
        {
          n =string1.charAt(i);
          if(isNaN(n)!= true){num=n;}
        }

        
        if(string1 == 'MFIRST0_value')
        {var type = sel.value;var e = document.getElementsByName("seq"+num)[0];
         
         e.style.display = 'none';
        
         var f = document.getElementsByName("cluster"+num)[0];
         
         f.style.display = 'none';var g = document.getElementsByName("pro"+num)[0];
         g.style.display = 'none';
         var h = document.getElementsByName("stringency")[0];
         h.style.display='none';
         }
          
        else
        {var type = sel[1].value;}
        
        if(string1 != 'MFIRST0_value')
        { var g=document.getElementsByName("pro"+num)[1];g.style.display = 'none';}
        
        box.options.length=0;
        box.style.display = 'block';
    
         var d=document.getElementsByName("pro"+num)[0];
        
         d.style.display = 'none';
         
         
        
        
        if (type == 'site_type')
        {%s}
        
        else if (type == 'GO_BP')
        {%s}
        
        else if (type == 'GO_MF')
        {%s}
        else if (type == 'GO_CC')
        {%s}
        
        else if (type == 'entrez_protein')
        {%s}
        else if (type == 'swissprot')
        {%s}
        else if (type == 'species')
        {%s}
        else if (type == 'pfam_site')
        {%s}
        else if (type == 'acc_gene')
        {
       
        %s
        box.style.display = 'none';
        d.disabled = false;
        d.style.display = 'block';
         g.disabled = false;
        g.style.display = 'block';
        
        }
        else if (type == 'pep_aligned')
        {
        d.disabled = false;
        d.style.display = 'block';
       
        box.style.display = 'none';
        }
        
        
        else if (type == 'gi')
        {%s}
        else if (type == 'domain')
        {%s}
        else if (type == 'scansite_bind')
        {
           %s
           h.style.display = 'block';
           h.disabled = false;
           

        }
        else if (type == 'scansite_kinase')
        {%s
        h.style.display = 'block';
           h.disabled = false;}
        else if (type == 'scansite')
        {%s}
        else if (type == 'clusters')
        {
           %s
           f.style.display = 'block';
           f.disabled = false;
           clusterBox(f,box);
           

        }
        else if (type=='pelm_kinase')
        {%s}
       
    }
    """ % (optionsDict['site_type'],optionsDict['GO_BP'],optionsDict['GO_MF'],optionsDict['GO_CC'],optionsDict['entrez_protein'],optionsDict['swissprot'],optionsDict['species'],optionsDict['pfam_site'],optionsDict['acc_gene'],optionsDict['gi'], optionsDict['domain'],optionsDict['scansite_bind'],optionsDict['scansite_kinase'],optionsDict['scansite'],optionsDict['clusters'],optionsDict['pelm'])
    print "</script>"
    print """<script type="text/javascript">"""
    if data:
        print """
   
    function addDataDiv(num,type)
    {
    
        if (type == "AND")
        {
            id = "AND"+num;
            var e = document.getElementById(id);
            var f = document.getElementById("OR"+num);
            e.style.display = 'block';
            document.getElementsByName(id)[0].value = "active";
            f.style.display = 'none';
            document.getElementsByName("OR"+num)[0].value = "inactive";
        }
        else
        {
            id = "OR"+num;
            var e = document.getElementById(id);
            var f = document.getElementById('AND'+num);
            e.style.display = 'block';
            f.style.display = 'none';
            document.getElementsByName(id)[0].value = "active";
            document.getElementsByName("AND"+num)[0].value = "inactive";
        }
        
    }
    function noDataDiv(num)
    {
        document.getElementById("AND"+num).style.display = 'none';
        document.getElementById("OR"+num).style.display = 'none';
        document.getElementsByName("AND"+num)[0].value = "inactive";
        document.getElementsByName("OR"+num)[0].value = "inactive";
    }"""
    print """
    function addMetaDataDiv(num,type)
    {
    
        if (type == "MAND")
        {
            id = "MAND"+num;
            var e = document.getElementById(id);
            var f = document.getElementById("MOR"+num);
            e.style.display = 'block';
            f.style.display = 'none';
            document.getElementsByName(id)[0].value = "active";
            document.getElementsByName("MOR"+num)[0].value = "inactive";
        }
        else
        {
            id = "MOR"+num;
            var e = document.getElementById(id);
            var f = document.getElementById("MAND"+num);
            e.style.display = 'block';
            f.style.display = 'none';
            document.getElementsByName(id)[0].value = "active";
            document.getElementsByName("MAND"+num)[0].value = "inactive";
        }
        
    }
    function noMetaDataDiv(num)
    {
        document.getElementById("MAND"+num).style.display = 'none';
        document.getElementById("MOR"+num).style.display = 'none';
        document.getElementsByName("MAND"+num)[0].value = "inactive";
        document.getElementsByName("MOR"+num)[0].value = "inactive";
    }
    """
    print """
    var active = 0;
    function filterDivs(form)
    {
         if(form.dataType.value == 'data')
         {
         document.getElementById('data').style.display = 'block';
         document.getElementById('metadata').style.display = 'none';
         }
         else if (form.dataType.value == 'metadata')
         {
         document.getElementById('metadata').style.display = 'block';
         document.getElementById('data').style.display = 'none';
         }
         else
         {
         document.getElementById('data').style.display = 'block';
         document.getElementById('metadata').style.display = 'block';
         }
    }</script>"""


    ### frames javascript
    print """<script type = 'text/javascript'>
    function changeTarget()
    {
        
        if(document.form1.select.value == 'all')
        {
        document.form1.target = "sub1";
        }

    }
    
    function reset()
    {
        document.getElementsByName('counter')[0].value = 0;
        document.form1.target = "sub1";
        for (i = 1; i <= 20; i = i+1)
        {
             document.getElementById("subset"+i).style.display = 'none';
             document.getElementsByName("b"+i)[0].src = "%stab_nonexistent.png";
             
        }
       
        document.form1.selectedids.value='';
        changeDivs(1);
    }
    function set()
    {
    //alert('in set');
      var win=document.getElementsByName('sub1')[0];
      //win.height=2000;
      document.form1.select.style.display = 'inline';
      var count = document.getElementsByName('counter')[0].value;
      if (count > 0){
      if (document.getElementsByName('sub'+count)[0].contentWindow.document.msform != null)
         {
         //alert('in set');
         var ids = document.getElementsByName('sub'+count)[0].contentWindow.document.msform.msids.value;
         document.getElementsByName('selectedids'+count)[0].value = ids;
         document.form1.selectedids.value = ids;
         //alert(document.form1.selectedids.value);
         
         }
      }
      if (document.form1.select.value != 'None')
      {
      if(document.form1.select.value == 'all'){document.form1.selectedids.value='';}
       var count = parseInt(document.getElementsByName('counter')[0].value);
       document.getElementsByName('counter')[0].value = count + 1;
       // set that frame to be active
       changeDivs(count + 1);
       newtarget = count+1;
     
       document.form1.target="sub"+newtarget;
       changeDivs(newtarget);
      
       
      
       
       }

       else
       {
       reset()
       }
      
    }
    
    function changeDivs(i)
    {
      
       for (j = 1; j <= %s; j = j+1)
       {
          
          var d = document.getElementById("subset" + i);
          level = document.getElementsByName("level" + j)[0].value;
         
          if (j == i)
          {
           
              
              document.getElementById("subset"+j).style.display = 'block';
              document.getElementsByName("level" + j)[0].value = 'active';
              
              document.getElementsByName("b"+j)[0].style.display = 'inline';
              
              document.getElementsByName("b"+j)[0].src = "%stab_active.jpg";
              document.getElementById("font"+j).style.display = 'block';
              document.getElementsByName('sub'+j)[0].style.height = document.getElementsByName('sub'+j)[0].contentWindow.document.height;
              //alert(document.getElementsByName('sub'+j)[0].contentWindow.document.style.height);
              
          }
          else if (level == 'nonexistent')
          {
              document.getElementById("subset"+j).style.display = 'none';
          }
          else 
          {
              
              document.getElementById("subset"+j).style.display = 'none';
              document.getElementsByName("level" + j)[0].value = "inactive";
              document.getElementsByName("b"+j)[0].src = "%stab_inactive.jpg";
              
              
          }
          newlevel = document.getElementsByName("level" + j)[0].value;
          
       }
     
     
       
    }
    function closeDiv(i)
    {
        document.getElementById('subset'+i).style.display = 'none';
        document.getElementsByName('level' + i)[0].value = 'nonexistent';
        document.getElementsByName('b'+i)[0].style.display = 'none';
        document.getElementById("font"+i).style.display = "none";
             
        var j = 1;
        var set = "False";
        
        for (j = i; j> 0; j--)
        {
        
        var activity = document.getElementsByName('level'+j)[0].value;
        if (activity != 'nonexistent')
        {
        changeDivs(j);
        set = "True";
        break;
        }
       
        
       }
     
                 
    }

    function clearForm()
    {
       var form = document.getElementsByName('form1')[0];
       form.reset();
       metaDataBox(document.getElementsByName('MFIRST0_value')[0],document.getElementsByName('MFIRST0')[0],form);
       for (i = 1; i < 7; i++){
       noDataDiv(i);
       noMetaDataDiv(i);
       }
    }
    
    
    """ % (imgPath,MAX_NUM_SUBSETS,imgPath,imgPath)
    
    print '</script>'
    ### end javascript functions


def getExperimentParameters(c,expid):
    ## get number of MSids
    
    query = "select * from MS where experiment_id=%s" % expid

    c.execute(query)
    x=c.fetchall()
    numMS = len(x)
    ## get number of proteins
    proteins = sets.Set([item[3] for item in x])
    numPro = len(proteins)
    ## get number of peptides
    query = "select * from MS_phosphopep join MS on MS_phosphopep.MS_id = MS.id where experiment_id=11;"
    c.execute(query)
    x=c.fetchall()
    numPeps = len(x)
    
    return numMS, numPro,numPeps
def cSort(stringList):
    """case-insensitive string comparison sort
    doesn't do locale-specific compare
    though that would be a nice addition
    usage: stringList = caseinsensitive_sort(stringList)"""

    tupleList = [(x.lower(), x) for x in stringList]
    tupleList.sort()
    return [x[1] for x in tupleList]


def printSelectOptionBox(list,GO=False,addNull=False):
    string =''
    start = 0
    if addNull:
        string +=  "box.options[0] = new Option('NULL','NULL');\n"
        start = 1
    if not(GO):
        list = cSort(list)
        for i in range(start,len(list)):
            if i < 10000:
                string += "box.options["+str(i) +"] = new Option('%s','%s');\n" % (list[i],list[i])
    else:
        list = sortListByElement(list,1)
        for i in range(start,len(list)):
            if i < 100000:
                string += "box.options["+str(i) +"] = new Option('%s','%s');" % (str(list[i][0])+' (' + str(list[i][1]).replace("'","")+')',list[i][0])
    if len(list)==0:
        string += "box.Options[" + str(0) + "] = new Option('Null','Null');"

    return string

def printClusterOptionBox(list):
    string = ''
    for item in list:
        string += 'box.options.add(new Option("Cluster Number: "+%s,%s));'% (str(item),str(item))
    return string
def getGO(exp_id,c):
    c.execute("""select protein_id from MS where experiment_id = %s""",(exp_id))
    x=c.fetchall()
    proteins=[item[0] for item in x]
    BP=[]
    MF=[]
    CC=[]
    for item in proteins:
        query = """select GO, aspect,term from GO join protein_GO on GO.id=protein_GO.GO_id where protein_GO.protein_id=%s sort by term""" % item
        c.execute(query)
        # make list of tuples (GO id, aspect)
        x=c.fetchall()
        GOs = [(term[0],term[1],term[2]) for term in x]
        for go in GOs:
            if go[1]=='P' and (go[0],go[2]) not in BP:
                BP.append((go[0],go[2]))
            elif go[1]=='F' and (go[0],go[2]) not in MF:
                MF.append((go[0],go[2]))
            else:
                if (go[0],go[2]) not in CC:append((go[0],go[2]))
    # uniquify lists
    #BP = reduce(lambda l, x: x not in l and l.append(x) or l, BP, [])
    #MF = reduce(lambda l, x: x not in l and l.append(x) or l, MF, [])
    #CC = reduce(lambda l, x: x not in l and l.append(x) or l, CC, [])
    return (BP,MF,CC)
## #make all of these take ms_id subset -> dict of (item: count)
## #can use for selection options and for enrichment, etc.

def getAcc(exp_id,c):
    c.execute("""select acc_gene from protein join MS on MS.protein_id=protein.id where experiment_id = %s""",(exp_id))
    x = c.fetchall()
    proteins= [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, proteins, [])

## returns a tuple of 3 lists (BP, MF, CC)
def getGO(exp_id,c):
    c.execute("""select protein_id from MS where experiment_id = %s""",(exp_id))
    x=c.fetchall()
    proteins=[item[0] for item in x]
    BP=[]
    MF=[]
    CC=[]
    for item in proteins:
        query = """select GO, aspect,term from GO join protein_GO on GO.id=protein_GO.GO_id where protein_GO.protein_id=%s""" % item
        c.execute(query)
        # make list of tuples (GO id, aspect)
        x=c.fetchall()
        GOs = [(term[0],term[1],term[2]) for term in x]
        for go in GOs:
            if go[1]=='P': BP.append((go[0],go[2]))
            elif go[1]=='F': MF.append((go[0],go[2]))
            else: CC.append((go[0],go[2]))
    # uniquify lists
    BP = reduce(lambda l, x: x not in l and l.append(x) or l, BP, [])
    MF = reduce(lambda l, x: x not in l and l.append(x) or l, MF, [])
    CC = reduce(lambda l, x: x not in l and l.append(x) or l, CC, [])
    return (BP,MF,CC)
#make all of these take ms_id subset -> dict of (item: count)
#can use for selection options and for enrichment, etc.

        
def getSpecies(exp_id,c):
    c.execute("""select species from protein join MS on MS.protein_id=protein.id where experiment_id=%s""",(exp_id))
    x= c.fetchall()
    species = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, species,[])

def getPfam(exp_id,c):
    c.execute("""select pfam_site from phosphopep join MS on MS.protein_id=phosphopep.protein_id where experiment_id=%s""",(exp_id))
    x= c.fetchall()
    pfam = [item[0] for item in x]
    pfam = reduce(lambda l, x: x not in l and l.append(x) or l, pfam,[])
    new = []
    for item in pfam:
        #if item != '~~~': new.append(item)
        new.append(item)
    return new

def getDomain(exp_id,c):
    c.execute("""select label from domain join MS on MS.protein_id= domain.protein_id where experiment_id=%s""",(exp_id,))
    x=c.fetchall()
    dom = [item[0] for item in x]
    dom = reduce(lambda l, x: x not in l and l.append(x) or l, dom,[])
    new = []
    for item in dom:
        #if item != '~~~': new.append(item)
        new.append(item)
    return new

def getEntrez(exp_id,c):
    c.execute("""select value from acc join MS on MS.protein_id=acc.protein_id where type = 'entrez_protein' and experiment_id = %s""",(exp_id))
    x=c.fetchall()
    entrez = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, entrez, [])

def getBind(exp_id,c):
    query = """select value from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id=MS_phosphopep.phosphopep_id join MS on MS.id=MS_phosphopep.MS_id where experiment_id = %s and source = 'scansite_bind'"""%exp_id
    c.execute(query)
    x=c.fetchall()
    binds = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, binds, [])
def getSKinase(exp_id,c):
    query = """select value from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id=MS_phosphopep.phosphopep_id join MS on MS.id=MS_phosphopep.MS_id where experiment_id = %s and source = 'scansite_kinase'"""%exp_id
    c.execute(query)
    x=c.fetchall()
    sks = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, sks, [])
def getScan(exp_id,c):
    query = """select value from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id=MS_phosphopep.phosphopep_id join MS on MS.id=MS_phosphopep.MS_id where experiment_id = %s and source = 'scansite' """%exp_id
    c.execute(query)
    x=c.fetchall()
    ks = [item[0] for item in x if item[0]!='~~~']
    return reduce(lambda l, x: x not in l and l.append(x) or l, ks, [])
def getPELM(exp_id,c):
    query = """select value from phosphopep_prediction join MS_phosphopep on phosphopep_prediction.phosphopep_id=MS_phosphopep.phosphopep_id join MS on MS.id=MS_phosphopep.MS_id where experiment_id = %s and source = 'pelm_kinase'"""%exp_id
    c.execute(query)
    x=c.fetchall()
    sks = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, sks, [])
def getGi(exp_id,c):
    c.execute("""select value from acc join MS on MS.protein_id=acc.protein_id where type = 'gi' and experiment_id = %s""",(exp_id,))
    x=c.fetchall()
    gi = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, gi, [])
def getSwiss(exp_id,c):
    c.execute("""select value from acc join MS on MS.protein_id=acc.protein_id where type = 'swissprot' and experiment_id = %s""",(exp_id))
    x=c.fetchall()
    swiss = [item[0] for item in x]
    return reduce(lambda l, x: x not in l and l.append(x) or l, swiss, [])


def printImageMap(i):
    print """
    <map name ="tabmap%s">
    <area shape = "rect" coords = "0,0,70,30" href="javascript:changeDivs(%s)" alt="make tab active">
    <area shape = "rect" coords = "70,0,100,30" href="javascript:closeDiv(%s)" alt="make tab inactive">
    </map>
    """ %(i,i,i)
def sortListByElement(list,index):
	tmpi = [item[index] for item in list]
        tmpi = cSort(tmpi)
        newlist = []
        for item in tmpi:
            for thing in list:
                if thing[index]==item:
                    newlist.append(thing)
        return newlist
	## tmp = [(item[index],item) for item in list]
	
## 	tmp.sort()
## 	list = [t[1] for t in tmp]
## 	return list

def getType(exp_id):
	# look at stuff for one ms_id, since all will (theoretically!) be the same
	c.execute("""select id from MS where experiment_id=%s""",(exp_id,))
	try:
		MS_id=c.fetchall()[0][0]
	except IndexError:
		return [],[]
	c.execute("""select * from data where MS_id=%s""",(MS_id,))
	x=c.fetchall()
	runtypes = [item[2] for item in x]
	# uniquify to get a list of all possible run types
	runtypes = reduce(lambda l, x: x not in l and l.append(x) or l, runtypes, [])
	# get types (same for all runs, so just find for one run type)
        try:
            run = runtypes[0]
        except IndexError:
            return [],[]
	# get types
	c.execute("""select type,label from data where MS_id=%s and run = %s order by priority""",(MS_id,run))
	x=c.fetchall()
	types = [(x[i][0],x[i][1]) for i in range(len(x))]
        if len(toAverage(runtypes))>1:
		runtypes.append('Averaged')
	return (runtypes, types)

main()

