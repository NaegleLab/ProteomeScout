from xml.dom.minidom import parse
import sys
import pickle

def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)

if __name__=='__main__':
    pfam_filename = 'scripts/data/PfamFamily.xml'

    sys.stderr.write( "Loading pfam XML definitions from %s\n" % (pfam_filename) )
    pfam_xml = parse(pfam_filename)

    db_nodes = pfam_xml.getElementsByTagName('database')
    sys.stderr.write( "Expecting %s entries\n" % ( getText(db_nodes[0].getElementsByTagName('entry_count')[0].childNodes )))

    pfam_mapping = {}

    for entry in db_nodes[0].getElementsByTagName('entries')[0].getElementsByTagName('entry'):
        name = getText( entry.getElementsByTagName('name')[0].childNodes )
        pfam_type = None

        addt_fields = entry.getElementsByTagName('additional_fields')[0]
        for field in addt_fields.getElementsByTagName('field'):
            if field.getAttribute('name') == 'type':
                pfam_type = getText( field.childNodes )

        if pfam_type != None:
            pfam_mapping[name] = pfam_type

    sys.stderr.write( "Got %d map entries\n" % ( len(pfam_mapping) ))

    print pickle.dumps(pfam_mapping)
