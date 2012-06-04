import template
def makeUrlDict():
    return {'swissprot':"http://ca.expasy.org/uniprot/",'entrez_protein':'http://www.ncbi.nlm.nih.gov/protein/','gi':'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=','refseq':'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id=','GO':'http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=','genbank':'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&amp;id='}

def makeGODict():
        return {'F':'FF00FF', 'C':'999900','B':'663399'}

def makeGOList(type, dictionary, url_dict, go_dict):
	print '<ul class="no_bullet">'
	for item in dictionary["go"]:		
		if item[0]==type:
			print "<li>"
			print """<a href='%s""" % url_dict.get('GO')
			print """%s' target='_blank'>""" % item[1]
			print """<font color= %s>""" % go_dict.get(item[0])
			print """ %s</font></a>""" % item[1]
			print """ (%s)""" % item[2]
			print "</li>"
	print '</ul>'


def makeProteinInfoDict(pid, c):
            d=dict()
            c.execute("""select * from protein where id=%s""",(pid,))
            x = c.fetchall()[0]
            sequence, species,acc_gene, name = x[1:5]
            c.execute("""select * from acc where protein_id=%s""",(pid,))
            x=c.fetchall()
            d["accs"]=[item for item in x] ## list of tuples of id, type, val, proid
            c.execute("""select aspect,GO,term from GO join protein_GO on GO.id=protein_GO.GO_id where protein_id=%s""",(pid,))
            x=c.fetchall()
            d["go"]=[item for item in x] ## list of tuples (id,GO_type,GO_id,name)
            d["sequence"]=sequence
            d["species"]=species
            d["acc"]=acc_gene
            d["name"]=name
            c.execute("""select * from domain where protein_id =%s""",(pid,))
            x=c.fetchall()
            d["domain_data"]=[item for item in x] ## list of tuples (id, label, start, stop, p_value, source, params, protein_id)
            return d


# Displays the detailed information about the protein
def makeProteinDetails(dictionary):
	url_dict = makeUrlDict()	
	print '<p><label>Name:</label>'
	print dictionary["name"]
	print '</p>'
	print '<p><label>Species:</label>'
	print dictionary["species"]
	print '</p>'
	print '<p><label>Accessions:</label>'
	print '<ul>'
	for item in dictionary["accs"]:
		print '<li><i>'
		print item[1] + '</i>: '
		if item[1]!='gene_synonym':
			print """<a href='%s""" % url_dict.get(item[1])
			print """%s' target='_blank'>""" % item[2]
		else:
			print str(item[2])
		if item[1]!='gene_synonym':
			print """%s</a>""" % item[2]
		print "</li>"
	
	print '</ul>'
	print '</p>'
	
def makeGeneOntologies(dictionary,c):
	url_dict = makeUrlDict()	
	go_dict = makeGODict()

	print '<div style="margin-bottom:1em">'
	template.printGODocs(c)
	print '</div>'
	
	print '<table cellspacing=0 cellpadding=0 id="gene_ontology_breakdown">'
	print '<tr>'
	print '<th>Molecular Functions</th>'
	print '<th>Cellular Components</th>'
	print '<th>Biological Process</th>'
	print '</tr>'
	
	
	print '<tr>'	
	print '<td valign="top">'
	makeGOList('F', dictionary, url_dict, go_dict)	
	print '</td>'

	print '<td valign="top">'
	makeGOList('C', dictionary, url_dict, go_dict)	
	print '</td>'

	print '<td valign="top">'
	makeGOList('P', dictionary, url_dict, go_dict)	
	print '</td>'
	print '</tr>'
	print '</table>'
	

	

def getDefault(pids,c):
    if len(pids)==0: return None
    """ return the default protein id from a list of protein ids by these criteria:
    1. Choose the isoform/protein with the most GO annotations
2. When there is a tie, choose the isoform/protein with the most assigned
modified peptides.
"""
    ## get number of GO annotations for each
    goAnnDict = {}
    numPepDict = {}
    for item in pids:
        pid = item[1][0]
        d = makeProteinInfoDict(pid, c)
        goAnnDict[pid]=len(d.get("go",[]))
        query = """select id from phosphopep where protein_id = %s"""%pid
        c.execute(query)
        x=c.fetchall()
        numPepDict[pid]=len(x)
    maxGO = max(goAnnDict.values())
    candidates = []
    for item in pids:
        pid = item[1][0]
        if goAnnDict[pid] == maxGO:
            candidates.append(pid)
    ## if only one, it's our default
    if len(candidates)==1:
        return candidates[0]
    maxPep = 0
    for item in candidates:
        if numPepDict[item] > maxPep: maxPep = numPepDict[item]
    for item in candidates:
        if numPepDict[item]==maxPep:
            return item
    return pids[0][1][0]
    
def getAllProteins(sequence,c,species):
    """ get all the possible ambiguous protein assignments for a given sequence and species"""
    c.execute("""select protein_id,name,acc_gene from ambiguity join protein on ambiguity.protein_id=protein.id where peptide=%s""",(sequence,))
    x = c.fetchall()
    ids = [(item[0],item[1],item[2]) for item in x]
    same_species_ids = []
    # make list of accs and remove duplicates
    accs=[id[2] for id in ids]
    accs = reduce(lambda l,x: x not in l and l.append(x) or l, accs, [])
    # get rid of protein ids for different species
    for thing in ids:
        id = thing[0]
        c.execute("""select species from protein where id = %s""",(id,))
        x=c.fetchall()[0][0]
        if x==species:
            same_species_ids.append(thing)
    # make list to return
    allProteins = []
    for item in accs:
        a = []
        for thing in same_species_ids:
            if thing[2]==item:
                a.append(thing)
        allProteins.append(a)
    return [item for item in allProteins if item !=[]]



        
