from ptmscout.database import DBSession, taxonomies
from DB_init import DatabaseInitialization
import traceback
from paste.deploy.loadwsgi import appconfig
import os, sys
import re



def parse_taxon_def(line):
    line = line.split("|")
    tax_id = int( line[0].strip() )
    parent_id = int( line[1].strip() )

    return tax_id, None if parent_id == tax_id else parent_id

if __name__=='__main__':
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'test.ini'))
        
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()


        infile = open('scripts/data/nodes.dmp', 'r')

        tax_id_map = {}

        sys.stderr.write( "Loading taxon hierarchy...\n" )
        for line in infile:
            tax_id, parent_id = parse_taxon_def(line)
            tax_id_map[tax_id] = parent_id


        sys.stderr.write( "Push taxons to DB...\n" )
        i = 0
        for tax_id in tax_id_map:
            node = taxonomies.getTaxonomyById(tax_id)
            parent_id = tax_id_map[tax_id]

            if node and parent_id:
                parent_node = taxonomies.getTaxonomyById(parent_id)
                while parent_node == None and parent_id != None:
                    parent_id = tax_id_map[parent_id]
                    parent_node = taxonomies.getTaxonomyById(parent_id)

                if parent_node:
                    node.parent_id = parent_node.node_id

                DBSession.add(node)

            i += 1
            if i % 10000 == 0:
                sys.stderr.write( "%d\n" % ( i ) )
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
