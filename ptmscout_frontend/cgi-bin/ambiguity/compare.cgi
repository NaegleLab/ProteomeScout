#!/usr/local/bin/python
from pathname import *
import cgi
import cgitb
if displayPythonErrors:
    cgitb.enable()
print "Content-Type: Text/html\n\n"
print """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">"""
import MySQLdb
import sys

sys.path.append(path+'includes')
sys.path.append(path+'imageScripts')
from template import *
from proteinFcns import *
from shapes2 import *
from graph import *
import checks
from proteinFcns import *
# Prints relative links to the various sections on this long page
# XXX BOGON ALERT - this is basically duplicate in browse/protein.cgi
def startSection(title, anchor, help=""):
    print """<h2><a name="%s">%s""" %(anchor,title)
    if help != "":
        print helper(help)
        
    print '</a></h2>'
    
    printNavigation(anchor)
    print '<div class="section">'

def endSection():
	print '</div>'

def printNavSection(anchor, title, display_location,last=False):
	if ( anchor == display_location ):
		print title
	else:
		print """<a href="#%s">%s</a>""" %(anchor, title)
	if (last != True ):
		print '|'
		

def printNavigation(display_location=""):
	print '<div class="protein_navigation">'
	printNavSection('details'	, 'Protein Information'		, display_location)
	printNavSection('expression'	, 'Expression'		, display_location)
	printNavSection('structure'	, 'Protein Structure'		, display_location, True)
	print "%s" %(helper("View_All_Ambiguity_For_a_Single_Peptide"))
	print '</div>'
# END BOGON ALERT


def main():
    form = cgi.FieldStorage()
    exp_id=form.getvalue("exp_id",None)
    if (exp_id == "None"):
        exp_id = None

    printHeader("Protein Ambiguity", exp_id)
    checks.printJS()

    db = MySQLdb.connect(user= "%s"%user, passwd="%s"%mysql, db="%s"%database)
    c=db.cursor()

    if form.has_key("peptide") and form["peptide"].value!="":
        sequence=form.getvalue("peptide")
    else:
        print "<BR>Error. No sequence found.<BR><BR>"
        sequence = ''
        return
        
    if form.has_key("species") and form["species"].value!='':
        species = form.getvalue("species")
    else:
        species = 'homo sapiens'
        
    
    allProteins = getAllProteins(sequence,c,species)

    print '<div id="content">'
    makeProteinInformation(allProteins,c,exp_id,species,form,sequence)
    makeExpression(allProteins,c,exp_id,species,form,sequence)
    makeProteinStructure(allProteins,c,exp_id,species,form,sequence)
    print '</div>'
    
    db.close()
    printFooter()

def makeProteinInformation(allProteins,c,exp_id,species,form,sequence):
    startSection("Protein Information", 'details')
    if (exp_id != None):
        print "<img height=20 src='%sstar.jpg' alt='star'> = Currently chosen protein<br><br>"%imgPath
            
    printTable(allProteins,c,exp_id,species,form,sequence)
    endSection()
            

def makeExpression(allProteins,c,exp_id,species,form,sequence):
    if form.has_key("probeid"):
        probes= form.getlist("probeid")
    else: probes = []

    # XXX BOGON ALERT - this seems duplicate in browse/protein.cgi
    table = form.getvalue("table","")
    columns = form.getvalue("columns","")
    ht,hc,mo = [],[],[]
    if form.has_key("human_tissue_cols"):
        ht = form.getlist("human_tissue_cols")
    if form.has_key("nci60_cols"):
        hc = form.getlist("nci60_cols")
    if form.has_key("mouse"):
        mo = form.getlist("mouse")
    if table=="tissue":
        columns = ht
    elif table=="cell":
        columns = hc
    elif table=="mouse":
        columns = mo
    else:
        columns = []
    # probe id stuff
    if table == "": table = "tissue"
    if columns == [] or 'all' in columns: columns = 'all'
    # END BOGON ALERT

    pids=[]
    for item in allProteins:
        for protein in item:
            pids.append(protein[0])
            # make list of lists of probe options
            # remove duplicates
        pids = reduce(lambda l, x: x not in l and l.append(x) or l, pids, [])
            
    probeOptions = []
    for pid in pids:
        if getProbeOptions(pid,c) != []:
            probeOptions.append(getProbeOptions(pid,c))
    defaultprobes = [item[0] for item in probeOptions]
    #uniquify probeids
    defaultprobes=reduce(lambda l, x: x not in l and l.append(x) or l, defaultprobes,[])

    startSection("Expression", "expression")
    

    
    if probes==[]:
        graphProbes = defaultprobes
    else:
        graphProbes = probes

    print """<BR><a target="_blank" href="%stitle=Expression_Data"><img border="0" alt="help" width="20" src="%shelp.jpg"></a>"""%(helpPath,imgPath)
    print "<h6>Expression Data comes from the Genomics Institute Of The Novartis Research Institute <a target='_blank' href='http://symatlas.gnf.org/'>SymAtlas Project</a></h6>"
    print "<h6><a href='%sdisclaimer.html' target='_blank'>GNF SymAtlas Copyright</a></h6>"%jsPath
    
    print '<a name = "graph"></a>'
    if len(graphProbes)>0:
        printForm(allProteins,sequence,species,c,exp_id,form)
            
        printSetForm(table,graphProbes,multiple=True)
        if table=='both':
            columnsBoth= [ht,hc]
            if ht==[]:
                columnsBoth[0]='all'
            if hc==[]:
                columnsBoth[1]='all'
                
            makeGraph(graphProbes,species,c,table='tissue',columns=columnsBoth[0],multiple=True)
            makeGraph(graphProbes,species,c,table='cell',columns=columnsBoth[1],multiple=True)
        else:
            makeGraph(graphProbes,species,c,table=table,columns=columns,multiple=True)
    else:
        print 'No probes Found<BR><BR><BR>'

    endSection()

def makeProteinStructure(allProteins,c,exp_id,species,form,sequence):
    startSection("Protein Structure", "structure", "Protein_Page#Protein_Structure")        
    
    printDomains(allProteins,c,exp_id,species,form,sequence)
    endSection()
        
    
                                                                                                                                 
## def getAllProteins(sequence,c,species):
##     c.execute("""select protein_id,name,acc_gene from ambiguity join protein on ambiguity.protein_id=protein.id where peptide=%s""",(sequence,))
##     x = c.fetchall()
##     ids = [(item[0],item[1],item[2]) for item in x]
##     same_species_ids = []
##     # make list of accs and remove duplicates
##     accs=[id[2] for id in ids]
##     accs = reduce(lambda l,x: x not in l and l.append(x) or l, accs, [])
##     # get rid of protein ids for different species
##     for thing in ids:
##         id = thing[0]
##         c.execute("""select species from protein where id = %s""",(id,))
##         x=c.fetchall()[0][0]
##         if x==species:
##             same_species_ids.append(thing)
##     # make list to return
##     allProteins = []
##     for item in accs:
##         a = []
##         for thing in same_species_ids:
##             if thing[2]==item:
##                 a.append(thing)
##         allProteins.append(a)
##     return [item for item in allProteins if item !=[]]

                                                                                    
def getChosen(pid,peptide,expid,c):
    query = """select experiment_id from MS where protein_id=%s and phosphopep = '%s'"""%(pid,peptide)
    c.execute(query)
    x=c.fetchall()
    if len(x)>0:
        exps = [str(item[0]) for item in x]
        if str(expid) in exps: return True
    return False


def printTable(allProteins, c,exp_id,species,form,sequence):

    url_dict = makeUrlDict()
    go_dict = makeGODict()
    
    print '<div style="margin-bottom:1em">'
    printGODocs(c)
    print '</div>'
    
    # start table of accessions
    print '<table id="protein_ambiguity" cellspacing=0 cellpadding=0>'
    # row one with contents (protein,accs, go)
    print '<tr><th>Protein</th><th>Accessions</th><th>GO Molecular Functions</th><th>GO Cellular Components</th><th>GO Biological Process</th></tr>'
    for item in allProteins:
        # print one row with the acc
        print '<tr bgcolor=99ccff><td colspan="5">%s</td></tr>' % item[0][2]
        # print other rows, 1 for each individual isoform
        for protein in item:
            ## if this is the protein they chose, put a star next to it
            chosen = getChosen(protein[0],sequence,exp_id,c)
            d = makeProteinInfoDict(protein[0],c)
            print '<tr valign=top >'
            print '<td>'
            
            if (chosen):
                print '<img height=20 src="%sstar.jpg" alt="chosen">' %(imgPath)                
            print '<a href = "%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">%s</a></td>' % (urlpath,str(protein[0]),str(exp_id),species,d["name"])

            print '<td>'
            
            # Print Accessions
            print '<ul class="no_bullet">'
            for item in d["accs"]:
                print '<li>'
                print item[1] + ': '
                if item[1]!='gene_synonym':
                    print """<a href='%s""" % url_dict.get(item[1])
                    print """%s'target='_blank'>""" % item[2]
                else: print '<a>'
                print """%s</a>""" % item[2]
                print "</li>"
            print '</td>'
            ## End Accessions

            molecularGODivID = 'F%s' % protein[0]
            cellularGODivId = 'C%s' % protein[0]
            biologicalGODivId = 'B%s' % protein[0]


            # print table of GO terms - functions table
            print '<td><div class="go_spacer">&nbsp;</div>'           
            print """<div id='%s' style="display:none">""" % molecularGODivID
            makeGOList('F', d, url_dict, go_dict)
            print '</div></td>'
            
            # print table of GO terms - cellular components
            print '<td><div class="go_spacer">'           
            
            print """<input type=button name="tog" value="Toggle GO terms" onClick = "toggleGO('%s','%s','%s',null)">""" % (molecularGODivID, cellularGODivId, biologicalGODivId)
            print '</div>'
            
            print '<div id="%s" style="display:none">' % cellularGODivId
            makeGOList('C', d, url_dict, go_dict)
            print '</div></td>'
             
            # print table of GO terms - biological process
            print '<td><div class="go_spacer">&nbsp;</div>'           
            print '<div id="%s" style="display:none">' % biologicalGODivId
            makeGOList('P', d, url_dict, go_dict)            
            print "</div></td>"
            # end table of GO terms
            
            print '</tr>'
    print '</table>'
    
def findMaxLength(proteins,c):
    #print proteins
    lengths=[]
    for p in proteins:
        lengths.append(getLength(p,c))
    if len(lengths)>0:
        return max(lengths)
    else:
        return 0
                      
def makePicture(pid,c,canvas,l,exp_id):
    d=makeDict(pid,c)
    p=makepDict(pid,c)
    olist,elist=makePOthersList(pid,c,exp_id)
    length = getLength(pid,c)
    p = Protein(pid,length,d,p,elist,olist,c,canvas=canvas,totallen=l)
    p.drawInJS(exp_id)
  

def printDomains(allProteins,c,exp_id,species,form,sequence):
    pids = []
    for item in allProteins:
        for thing in item:
            pids.append(thing[0])
    maxlen = findMaxLength(pids,c)
    for acc in allProteins:
        print "<BR><BR><h2>%s</h2>"%acc[0][2]
        
        for protein in acc:
            pid=protein[0]
            #print pid
            olist,elist=makePOthersList(pid,c,exp_id)
            
            chosen = getChosen(protein[0],sequence,exp_id,c)
            if not(chosen):
                print '<h5 align=center><a href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">%s</a></h5>'% (urlpath,str(protein[0]),str(exp_id),species,protein[1][:-1])
            else:
                print '<h5 align=center><img src="%sstar.jpg" alt="chosen" height=20><a href="%sbrowse/protein.cgi?protein_id=%s&amp;exp_id=%s&amp;species=%s">%s</a></h5>'% (imgPath,urlpath,str(protein[0]),str(exp_id),species,protein[1][:-1])
            print '('+getAcc(protein[0],c)+')'+'<br><BR><BR>'
            if len(olist)==0:
                print """<center><div id="%s" style='position:relative;left:200px;' ></center>"""%protein[1]
            else:
                print """<center><div id="%s" style='position:relative;left:200px;'></center><BR>"""%protein[1]
            pid = protein[0]
            pidlen = getLength(pid,c)
            totallen = pidlen*1.0/maxlen*900
            makePicture(protein[0],c,protein[1],totallen,exp_id)
            print "</div><BR><BR><BR>"
            print "<br>"*20

def printForm(allProteins,peptide,species,c,exp_id,form):
    pids = []
    names = []
    for item in allProteins:
        for protein in item:
            pids.append(protein[0])
            names.append((protein[2],getAcc(protein[0],c)))
    data='tissue'
    allProbes = []
    for protein in pids:
        if getChosen(protein,peptide,exp_id,c):
            allProbes.append(getProbeOptions(protein,c)+["chosen"])
        else:
            allProbes.append(getProbeOptions(protein,c))
    # make list of tuples (probeset list, [protein names])
    groupSameProbes = []
    for item in allProbes:
        if item not in [thing[0] for thing in groupSameProbes]:
            proName = []
            for i in range(len(allProbes)):
                if allProbes[i] == item: proName.append(names[i])
            groupSameProbes.append((item, proName))
    
    print '<center><BR><BR><BR>'
    width = min(len(groupSameProbes)*50,80)+10
    print '<div style="background-color:#ccccff;border:1px gray solid;width:%s%%;padding:5px;">'%width
        # fill in action
    print '<form name = "form" method=POST action="%sambiguity/compare.cgi#exp_display">'%(urlpath)
    
    
    
     
 
    print '<p align=left>'
    # table with grouping same probe sets
    print "Note: all gene products listed in a single column share a probeid"
    print '<p align=left>'
    print "Accessions listed are those used in ambiguity resolution<BR>"
    print '<p align=left>'
    print "GP=Gene Product"
    
    print '<table cellpadding = "5" ><tr>'
    if groupSameProbes !=[]:
        maxnum = max([len(item[1]) for item in groupSameProbes])
    else:
        maxnum=0
        print 'No probes available'
    usedaccs=[]
    countaccs=[]
    for i in range(maxnum):
        print '<tr>'
        for p in groupSameProbes:
            
            print '<td>'
            if len(p[1]) > i:
                try:
                    if p[0][-1]=="chosen":
                        chosen = True
                        p[0].remove("chosen")
                    else: chosen = False
                
                except IndexError: chosen = False
                countaccs.append(p[1][i])
                if p[1][i] not in usedaccs:
                    if chosen: print "<img height=20 src='%sstar_blue.jpg' alt='star'>"%imgPath
                    print p[1][i][0]
                    
                    
                    print "GP %s"%(countaccs.count(p[1][i]))
                    print "<BR>"
                    print "("+p[1][i][1] +")"
                    usedaccs.append(p[1][i])
                    
                else:
                    print p[1][i]+" Gene Product %s"%(countaccs.count(p[1][i]))
                    
            print '</td>'
        print '</tr>'
    for i in range(len(groupSameProbes)):
        print '<td>'
        if len(groupSameProbes[i][0])>0:
            print '<select name = "probeid">'
            for probe in groupSameProbes[i][0]:
                 print '<option value = "%s">%s</option>' %(probe,probe)
            print '</select>'
        else:
            print 'No probes found'
        print '</td>'   
    print '</tr></table>'
    print '<br><br>'
    print 'Data Type:'
    print '<select name = "table" onChange="changeChecks(this.form)">'
    if species =="homo sapiens":
        print '<option value="cell">NCI60</option>'
        print '<option value="tissue">Human Expression</option>'
        print '<option value="both">View Both</option>'
        if data =='tissue': cols=EXPRESSION_HUMAN_COLS
        else: cols=EXPRESSION_NCI60_COLS
    else:
        print '<option value = "null">Mouse Expression</option>'
        cols = EXPRESSION_MOUSE_COLS
    print '</select>'
    print '<p align=left><br>'
    print '<input type = hidden name="species" value = "%s">'% species
    print 'Select columns to display:'
    print '<input type=button name="buttons" value="show column choices" onClick="addRemoveChecks(this.form)"><br>'
    print '<input type="hidden" name="peptide" value=%s>' % peptide
    print '<br><INPUT TYPE="submit" NAME="button" Value="Display">'
    # begin checkbox div

    print '<a name = "exp_display"></a>'

    print '<div id="checks">'
    print '</div>'
    print '<input type="hidden" name="exp_id" value = %s>' % str(exp_id)
    print '</form>'
    print '</div>'
    print '</center>'
                                             
     
                        
    

main()
