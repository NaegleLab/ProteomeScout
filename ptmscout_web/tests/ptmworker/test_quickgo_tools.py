from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker import quickgo_tools

class IntegrationTestQuickGoQuery(IntegrationTestCase):
    
    def test_OBO_XML_Parser(self):
        oboxml = \
"""
<obo>
    <header>
        <format-version>1.2</format-version>
        <auto-generated-by>QuickGO: http://www.ebi.ac.uk/QuickGO</auto-generated-by>
        <synonymtypedef>
            <id>systematic_synonym</id>
            <name>Systematic synonym</name>
            <scope>EXACT</scope>
        </synonymtypedef>
        <default-namespace>gene_ontology</default-namespace>
        <remark>
            This is not the master copy of the OBO files, see http://www.geneontology.org/GO.downloads.ontology.shtml
        </remark>
    </header>
    <term>
        <id>GO:0005768</id>
        <name>endosome</name>
        <namespace>cellular_component</namespace>
        
        <def>
            <defstr>
            A membrane-bounded organelle to which materials ingested by endocytosis are delivered.
            </defstr>
        </def>
        
        <xref>
            <acc>sao1720343330</acc>
            <dbname>NIF_Subcellular</dbname>
        </xref>
        <xref>
            <acc>Endosome</acc>
            <dbname>UniProtKB-KW</dbname>
        </xref>
        <xref>
            <acc>Endosome</acc>
            <dbname>UniProtKB-SubCell</dbname>
        </xref>
        <xref>
            <acc>Endosome</acc>
            <dbname>Wikipedia</dbname>
        </xref>
        
        <is_a>GO:0044444</is_a>
        <is_a>GO:0043231</is_a>
    </term>
</obo>"""

        oboxml = quickgo_tools.OBOXMLParser(oboxml)
        
        self.assertEqual(1, len(oboxml.entries))
        obo = oboxml.entries[0]
        
        self.assertEqual("1.2", oboxml.version)
        self.assertEqual("GO:0005768", obo.goId)
        self.assertEqual("endosome", obo.goName)
        self.assertEqual('C', obo.goFunction)
        self.assertEqual(["GO:0044444","GO:0043231"], obo.is_a)
        

    def test_get_GO_term_should_return_term_data(self):
        version, obo = quickgo_tools.get_GO_term("GO:0005768")
        
        self.assertEqual("1.2", version)
        self.assertEqual("GO:0005768", obo.goId)
        self.assertEqual("endosome", obo.goName)
        self.assertEqual('C', obo.goFunction)
        self.assertEqual(["GO:0044444","GO:0043231"], obo.is_a)
    
    def test_batch_get_GO_annotations_should_return_annotation_subset(self):
        annotations, gene_symbols = quickgo_tools.batch_get_GO_annotations(['P50914', 'Q07912', 'Q8N9T8'])
        
        self.assertFalse('GO:0005622' in annotations['P50914'])
        self.assertEqual('20121125', annotations['P50914']['GO:0005515'])
        self.assertEqual('20121125', annotations['Q07912']['GO:0005515'])
        self.assertEqual({}, annotations['Q8N9T8'])
        

        self.assertEqual('TNK2', gene_symbols['Q07912'])
        self.assertEqual('RPL14', gene_symbols['P50914'])
