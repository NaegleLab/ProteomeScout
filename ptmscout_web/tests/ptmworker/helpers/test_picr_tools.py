from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.helpers import picr_tools

class IntegrationTestPICR(IntegrationTestCase):
    
    def test_picr_parser(self):
        picrxml = \
"""
<ns2:getUPIForAccessionResponse xmlns="http://model.picr.ebi.ac.uk" xmlns:ns2="http://www.ebi.ac.uk/picr/AccessionMappingService">
    <ns2:getUPIForAccessionReturn containsActiveSPTRCrossReference="false">
        <CRC64>F308F907C2F899E0</CRC64>
        <UPI>UPI0001838837</UPI>
        <identicalCrossReferences>
            <accession>IPI01014284</accession>
            <accessionVersion>1</accessionVersion>
            <databaseDescription>International Protein Index</databaseDescription>
            <databaseName>IPI</databaseName>
            <dateAdded xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <dateDeleted xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <deleted>false</deleted>
            <gi xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <taxonId xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        </identicalCrossReferences>
        <sequence xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        <timestamp>2009-01-12T00:00:00.000Z</timestamp>
    </ns2:getUPIForAccessionReturn>
    <ns2:getUPIForAccessionReturn containsActiveSPTRCrossReference="true">
        <CRC64>FCF6DC22DEF001DF</CRC64>
        <UPI>UPI0000DB2E73</UPI>
        <identicalCrossReferences>
            <accession>P29375</accession>
            <accessionVersion>3</accessionVersion>
            <databaseDescription>UniProtKB/Swiss-Prot</databaseDescription>
            <databaseName>SWISSPROT</databaseName>
            <dateAdded xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <dateDeleted xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <deleted>false</deleted>
            <gi xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <taxonId xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        </identicalCrossReferences>
        <logicalCrossReferences>
            <accession>IPI01014284</accession>
            <accessionVersion xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <databaseDescription>IPI</databaseDescription>
            <databaseName>IPI</databaseName>
            <dateAdded xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <dateDeleted xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <deleted>false</deleted>
            <gi xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <taxonId xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        </logicalCrossReferences>
        <logicalCrossReferences>
            <accession>IPI00021363</accession>
            <accessionVersion xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <databaseDescription>IPI</databaseDescription>
            <databaseName>IPI</databaseName>
            <dateAdded xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <dateDeleted xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <deleted>false</deleted>
            <gi xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <taxonId xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        </logicalCrossReferences>
        <sequence xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        <timestamp>2009-01-12T00:00:00.000Z</timestamp>
    </ns2:getUPIForAccessionReturn>
    <ns2:getUPIForAccessionReturn containsActiveSPTRCrossReference="false">
        <CRC64>8CFF8A88AE69A652</CRC64>
        <UPI>UPI000013318B</UPI>
        <identicalCrossReferences>
            <accession>IPI00021363</accession>
            <accessionVersion>3</accessionVersion>
            <databaseDescription>International Protein Index</databaseDescription>
            <databaseName>IPI</databaseName>
            <dateAdded xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <dateDeleted xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <deleted>false</deleted>
            <gi xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
            <taxonId xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        </identicalCrossReferences>
        <sequence xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true"/>
        <timestamp>2009-01-12T00:00:00.000Z</timestamp>
    </ns2:getUPIForAccessionReturn>
</ns2:getUPIForAccessionResponse>
"""
        response = picr_tools.PICRParser(picrxml)
            
        self.assertEqual(6, len(response.references))
        self.assertEqual(('ipi', 'IPI01014284'), response.references[0])
        self.assertEqual(('ipi', 'IPI01014284.1'), response.references[1])
        self.assertEqual(('swissprot', 'P29375'), response.references[2])
        self.assertEqual(('swissprot', 'P29375.3'), response.references[3])
        self.assertEqual(('ipi', 'IPI00021363'), response.references[4])
        self.assertEqual(('ipi', 'IPI00021363.3'), response.references[5])

    def test_get_picr_3(self):
        r1 = picr_tools.get_picr('gi|86168', 9031)
        r2 = picr_tools.get_picr('P03891', 9606)

        self.assertEqual(set([('ipi', 'IPI00596848'), ('ipi', 'IPI00596848.1'), ('swissprot', 'P68034'), ('swissprot', 'P68034.1')]), set(r1))
        self.assertEqual(set([('swissprot', 'P03891'), ('swissprot', 'P03891.2'), ('refseq', 'YP_003024027'), ('refseq', 'YP_003024027.1'), ('ipi', 'IPI00007979'), ('ipi', 'IPI00007979.4')]), set(r2))

    def test_get_picr_2(self):
        result = picr_tools.get_picr('EAX03173', 9606)

        self.assertEqual([], result)

    def test_get_picr(self):
        result = picr_tools.get_picr('P29375', 9606)

        self.assertEqual(8, len(result))
        self.assertTrue(('swissprot', 'P29375') in result)
        self.assertTrue(('swissprot', 'P29375.3') in result)
        self.assertTrue(('refseq', 'NP_001036068') in result)
        self.assertTrue(('refseq', 'NP_001036068.1') in result)
        self.assertTrue(('ipi', 'IPI01014284') in result)
        self.assertTrue(('ipi', 'IPI01014284.1') in result)
        self.assertTrue(('ipi', 'IPI00021363') in result)
        self.assertTrue(('ipi', 'IPI00021363.3') in result)
        
