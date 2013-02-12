from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.helpers import scansite_tools

class IntegrationTestScansiteQuery(IntegrationTestCase):
    
    def test_scansite_parser(self):
        scansitexml = \
"""<proteinScanResult>
    <predictedSite>
        <motifName>p85 SH2</motifName>
        <motifNickName>p85_SH2</motifNickName>
        <score>0.6751</score>
        <site>Y8</site>
        <siteSequence>LKKVVALyDYMPMNA</siteSequence>
    </predictedSite>
    <predictedSite>
        <motifName>Shc SH2</motifName>
        <motifNickName>Shc_SH2</motifNickName>
        <score>0.5591</score>
        <site>Y8</site>
        <siteSequence>LKKVVALyDYMPMNA</siteSequence>
    </predictedSite>
    <predictedSite>
        <motifName>Insulin Receptor Kinase</motifName>
        <motifNickName>InsR_Kin</motifNickName>
        <score>0.6306</score>
        <site>Y10</site>
        <siteSequence>KVVALYDyMPMNA**</siteSequence>
    </predictedSite>
    <proteinName>PTMSCOUT_QUERY</proteinName>
    <proteinSequence>LKKVVALYDYMPMNA</proteinSequence>
</proteinScanResult>"""
        parser = scansite_tools.ScansiteParser(scansitexml)
        
        self.assertEqual(3, len(parser.sites))
        
        self.assertEqual("p85 SH2", parser.sites[0].name)
        self.assertEqual("p85_SH2", parser.sites[0].nickname)
        self.assertEqual(0.6751, parser.sites[0].score)
        self.assertEqual("LKKVVALyDYMPMNA", parser.sites[0].sequence)
        
        self.assertEqual("Shc SH2", parser.sites[1].name)
        self.assertEqual("Shc_SH2", parser.sites[1].nickname)
        self.assertEqual(0.5591, parser.sites[1].score)
        self.assertEqual("LKKVVALyDYMPMNA", parser.sites[1].sequence)
        
        self.assertEqual("Insulin Receptor Kinase", parser.sites[2].name)
        self.assertEqual("InsR_Kin", parser.sites[2].nickname)
        self.assertEqual(0.6306, parser.sites[2].score)
        self.assertEqual("KVVALYDyMPMNA**", parser.sites[2].sequence)
       
    def test_get_scansite_motif_with_whole_protein_should_get_everythin(self):
        prot = 'MAAVILESIFLKRSQQKKKTSPLNFKKRLFLLTVHKLSYYEYDFERGRRGSKKGSIDVEKITCVETVVPEKNPPPERQIPRRGEESSEMEQISIIERFPYPFQVVYDEGPLYVFSPTEELRKRWIHQLKNVIRYNSDLVQKYHPCFWIDGQYLCCSQTAKNAMGCQILENRNGSLKPGSSHRKTKKPLPPTPEEDQILKKPLPPEPAAAPVSTSELKKVVALYDYMPMNANDLQLRKGDEYFILEESNLPWWRARDKNGQEGYIPSNYVTEAEDSIEMYEWYSKHMTRSQAEQLLKQEGKEGGFIVRDSSKAGKYTVSVFAKSTGDPQGVIRHYVVCSTPQSQYYLAEKHLFSTIPELINYHQHNSAGLISRLKYPVSQQNKNAPSTAGLGYGSWEIDPKDLTFLKELGTGQFGVVKYGKWRGQYDVAIKMIKEGSMSEDEFIEEAKVMMNLSHEKLVQLYGVCTKQRPIFIITEYMANGCLLNYLREMRHRFQTQQLLEMCKDVCEAMEYLESKQFLHRDLAARNCLVNDQGVVKVSDFGLSRYVLDDEYTSSVGSKFPVRWSPPEVLMYSKFSSKSDIWAFGVLMWEIYSLGKMPYERFTNSETAEHIAQGLRLYRPHLASEKVYTIMYSCWHEKADERPTFKILLSNILDVMDEES'
        sites = scansite_tools.get_scansite_motif(prot, "MAMMALIAN", filter_exact=False)

        self.assertEqual(134, len(sites))

    def test_get_scansite_motif_should_return_relevant_motif_data_for_partial_peps(self):
        pep = 'KVVALYDyMPMNA  '
        sites = scansite_tools.get_scansite_motif(pep, "MAMMALIAN")

        self.assertEqual(0, len(sites))

    def test_get_scansite_motif_should_return_relevant_motif_data(self):
        sites = scansite_tools.get_scansite_motif("LKKVVALyDYMPMNA", "MAMMALIAN")
                
        self.assertEqual(2, len(sites))
        
        self.assertEqual("p85 SH2", sites[0].name)
        self.assertEqual("p85_SH2", sites[0].nickname)
        self.assertEqual(0.6751, sites[0].score)
        self.assertEqual("LKKVVALyDYMPMNA", sites[0].sequence)
        
        self.assertEqual("Shc SH2", sites[1].name)
        self.assertEqual("Shc_SH2", sites[1].nickname)
        self.assertEqual(0.5591, sites[1].score)
        self.assertEqual("LKKVVALyDYMPMNA", sites[1].sequence)
