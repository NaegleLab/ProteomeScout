from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker import picr_tools

class IntegrationTestScansiteQuery(IntegrationTestCase):
    
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
            
        self.assertEqual(3, len(response.references))
        self.assertEqual(('ipi', 'IPI01014284', '1'), response.references[0])
        self.assertEqual(('swissprot', 'P29375', '3'), response.references[1])
        self.assertEqual(('ipi', 'IPI00021363', '3'), response.references[2])
    
    def test_get_picr_2(self):
        result = picr_tools.get_picr('EAX03173', 9606)

        self.assertEqual([], result)

    def test_get_picr(self):
        result = picr_tools.get_picr('P29375', 9606)

        self.assertEqual(4, len(result))
        self.assertEqual(('swissprot', 'P29375', '3'), result[0])
        self.assertEqual(('refseq', 'NP_001036068', '1'), result[1])
        self.assertEqual(('ipi', 'IPI01014284', '1'), result[2])
        self.assertEqual(('ipi', 'IPI00021363', '3'), result[3])
        
