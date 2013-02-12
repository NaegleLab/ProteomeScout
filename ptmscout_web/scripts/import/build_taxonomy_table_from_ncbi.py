from ptmscout.database import DBSession, taxonomies
from DB_init import DatabaseInitialization
import traceback
import os, sys
import re


def parse_node_types(line):
    line = line.split("|")
    tax_id = int( line[0].strip() )
    node_type = line[2].strip()

    return tax_id, node_type

def parse_taxon_def(line):
    line = line.split("|")
    tax_id = int( line[0].strip() )
    name = line[1].strip()
    tp = line[3].strip()

    return tax_id, name, tp


def get_scientific_name(name_map):
    for n in name_map:
        if name_map[n] == 'scientific name':
            return n


if __name__=='__main__':
    try:
        settings = os.path.join('data', 'ptmscout', 'ptmscout_web', 'production.ini')

        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        sys.stderr.write( "Loading node types...\n" )
        node_file = open('scripts/data/nodes.dmp', 'r')

        skip_types = set(['superkingdom', 'species', 'no rank', 'kingdom'])
        tax_type_map = {}
        for line in node_file:
            tax_id, node_type = parse_node_types(line)
            tax_type_map[tax_id] = node_type


        sys.stderr.write( "Loading taxon node names...\n" )
        infile = open('scripts/data/names.dmp', 'r')
        tax_id_map = {}
        for line in infile:
            tax_id, name, tp = parse_taxon_def(line)
            name_map = tax_id_map.get(tax_id, {})
            name_map[name] = tp
            tax_id_map[tax_id] = name_map


        sys.stderr.write( "Push taxons to DB...\n" )
        i = 0
        created = 0
        for tax_id in tax_id_map:
            if tax_type_map[tax_id] not in skip_types:
                node = taxonomies.getTaxonomyById(tax_id)

                if not node:
                    scientific_name = get_scientific_name(tax_id_map[tax_id])

                    if scientific_name != None:
                        node = taxonomies.Taxonomy()
                        node.node_id = tax_id
                        node.kingdom = ''
                        node.name = scientific_name
                        node.strain = None
                        node.parent_id = None
                        created += 1
                        DBSession.add(node)

            i += 1
            if i % 10000 == 0:
                sys.stderr.write( "%d, Created: %d\n" % ( i, created ) )
                created = 0
                DBSession.flush()
                dbinit.commit()
                dbinit.new_transaction()


    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
