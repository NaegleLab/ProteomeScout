from tests.PTMScoutTestCase import IntegrationTestCase
from ptmworker.helpers import pfam_tools
from ptmscout.database import protein
from mock import patch

class IntegrationTestPFamQuery(IntegrationTestCase):
    
    def test_pfam_parser(self):
        pfam_xml = \
"""<?xml version="1.0" encoding="UTF-8"?>
<pfam xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xmlns="http://pfam.sanger.ac.uk/"
      xsi:schemaLocation="http://pfam.sanger.ac.uk/
                          http://pfam.sanger.ac.uk/static/documents/schemas/results.xsd"
      release="26.0" 
      release_date="2011-11-17">
  <results job_id="8b00a2c5-ccef-4e55-b7cf-c1f9c500b21a">
    <matches>
      <protein length="280">
        <database id="pfam" release="26.0" release_date="2011-11-17">
          <match accession="PF13465.1" id="zf-H2C2_2" type="Pfam-A" class="Domain">
            <location start="185" end="210" ali_start="188" ali_end="209" hmm_start="4" hmm_end="25" evalue="0.00093" bitscore="19.2" evidence="hmmer v3.0" significant="1">
              <hmm>
                <![CDATA[rHlrkHtgekpyeCplCgktFk]]>
              </hmm>
              <match_string>
                <![CDATA[ H+r+H+g+k  eC  C+kt+ ]]>
              </match_string>
              <pp>
                <![CDATA[7******************975]]>
              </pp>
              <seq>
                <![CDATA[QHQRIHSGKKLCECNECSKTCN]]>
              </seq>
              <raw>
                <![CDATA[
#HMM       rHlrkHtgekpyeCplCgktFk
#MATCH      H+r+H+g+k  eC  C+kt+ 
#PP        7******************975
#SEQ       QHQRIHSGKKLCECNECSKTCN
                ]]>
              </raw>
            </location>
            <location start="213" end="237" ali_start="215" ali_end="234" hmm_start="3" hmm_end="22" evalue="0.00094" bitscore="19.2" evidence="hmmer v3.0" significant="1">
              <hmm>
                <![CDATA[rrHlrkHtgekpyeCplCgk]]>
              </hmm>
              <match_string>
                <![CDATA[ rH r+H+gekpyeC++ gk]]>
              </match_string>
              <pp>
                <![CDATA[59*************97555]]>
              </pp>
              <seq>
                <![CDATA[IRHHRIHSGEKPYECDKHGK]]>
              </seq>
              <raw>
                <![CDATA[
#HMM       rrHlrkHtgekpyeCplCgk
#MATCH      rH r+H+gekpyeC++ gk
#PP        59*************97555
#SEQ       IRHHRIHSGEKPYECDKHGK
                ]]>
              </raw>
            </location>
            <location start="241" end="266" ali_start="241" ali_end="264" hmm_start="1" hmm_end="24" evalue="0.017" bitscore="15.2" evidence="hmmer v3.0" significant="1">
              <hmm>
                <![CDATA[nLrrHlrkHtgekpyeCplCgktF]]>
              </hmm>
              <match_string>
                <![CDATA[+L+rHl++H g k +eC+  g+ F]]>
              </match_string>
              <pp>
                <![CDATA[599**************9766555]]>
              </pp>
              <seq>
                <![CDATA[DLVRHLKIHIGDKVLECSERGRIF]]>
              </seq>
              <raw>
                <![CDATA[
#HMM       nLrrHlrkHtgekpyeCplCgktF
#MATCH     +L+rHl++H g k +eC+  g+ F
#PP        599**************9766555
#SEQ       DLVRHLKIHIGDKVLECSERGRIF
                ]]>
              </raw>
            </location>
          </match>
          <match accession="PF13900.1" id="GVQW" type="Pfam-A" class="Domain">
            <location start="61" end="112" ali_start="66" ali_end="112" hmm_start="7" hmm_end="51" evalue="0.024" bitscore="14.7" evidence="hmmer v3.0" significant="0">
              <hmm>
                <![CDATA[AGVqWhdLgSLQ..pppPrfkr..fsclslpsswdyrhhpprlanFvff]]>
              </hmm>
              <match_string>
                <![CDATA[A+ ++    S +  +++P+++   +++++l + +++ hh+   ++Fvf+]]>
              </match_string>
              <pp>
                <![CDATA[6666666667666667788877558999***********99..****97]]>
              </pp>
              <seq>
                <![CDATA[ARLECSGPISAHcnLHLPGSSNspAPASQLAGITGTHHHTR--LIFVFL]]>
              </seq>
              <raw>
                <![CDATA[
#HMM       AGVqWhdLgSLQ..pppPrfkr..fsclslpsswdyrhhpprlanFvff
#MATCH     A+ ++    S +  +++P+++   +++++l + +++ hh+   ++Fvf+
#PP        6666666667666667788877558999***********99..****97
#SEQ       ARLECSGPISAHcnLHLPGSSNspAPASQLAGITGTHHHTR--LIFVFL
                ]]>
              </raw>
            </location>
          </match>
        </database>
      </protein>
    </matches>
  </results>
</pfam>
"""
        pfam_parser = pfam_tools.PFamParser(pfam_xml)
        self.assertEqual(4, len(pfam_parser.domains))
        
        testdom = pfam_parser.domains[0]
        self.assertEqual(185, testdom.start)
        self.assertEqual(210, testdom.stop)
        self.assertEqual(1, testdom.significant)
        self.assertEqual("zf-H2C2_2", testdom.label)
        self.assertEqual("PF13465.1", testdom.accession)
        self.assertEqual(0.00093, testdom.p_value)
        self.assertEqual("26.0", testdom.release)
        self.assertEqual("Domain", testdom.class_)
        
        testdom = pfam_parser.domains[3]
        self.assertEqual(61, testdom.start)
        self.assertEqual(112, testdom.stop)
        self.assertEqual(0, testdom.significant)
        self.assertEqual("GVQW", testdom.label)
        self.assertEqual("PF13900.1", testdom.accession)
        self.assertEqual(0.024, testdom.p_value)
        self.assertEqual("26.0", testdom.release)
        self.assertEqual("Domain", testdom.class_)
        
    
    def test_get_computed_pfam_domains_3(self):
        prot_seq = "MASWESRKLLLLLWKNFTLKRRKFGTLVSEIVLVLLLSIVLLTTRHLLSIKKIEALYFPDQPISTVPSFFRASVTSFYPWELAYVPSNIMVVENIVKNVKNDLNLHMKVIGFPSESDFEDYARSTVNSRNILAAIVFGHNFANSSDPLPKKVKYYLRFSDIKKNINSGAYYQGDTWLTKFLFHSLRLVGPRNPYEADGGSPGYITEGFLAVQHALDKAIMLHHGGADAAALFNDISLFIQRFPYPAYYHDYFYLFATTFIPLTVACTFFFNHYVLVWSIVWEKENRLKEYQLMIGLRNWMFWVAYFFTFLCLYFINIIVMCMVLFVKIEPAPIFQYNDPTLVFIFLLFYAISSIFFSFMVSTLFNKVSLAMSLGSFLFFLTYFPAVAMHQSFERMPSKQKLIWSFDFNVGMAFGFRFLVNTDAKKTGMKWSNIFLSTDSDSFLFAYVLGMLLADAFIYGLVAWYIEAVFPGEYGVPKPWNFFLMHSYWFGEPPQQKLEITQFYERVESKYFEAEPTDLTAGIQIKHLHKVFQKNNTTKVAIKDLSLNLYEGQVTVLLGHNGAGKSTTLSILSGLYPPTSGEAYVHGEDISQHMDQVRNSLGLCPQQNLLFDHLTVSEHLYFYCRIKGVPQKMYLEETNNMLSAFNLMEKCDAFSKSLSGGMKRKLAIIIALIGGSKVAILDEPTSGMDPASRRSTWDILQTYKQNRTILLTTHYMDEADVLGDRIAIMVRGTLRCCGSSVFLKRLYGVGSHLVMVKEPYCDIAEISKLIHSYVPTATLETNVGNELSFILPKEYTHRFEALFTALEENQENLGISSFGVSITTMEEVFLKVSNLEDSKTDIEATQSPSVGSKGNKNGDVESSGRVGFPTQSEDQNIVFNTGCSLYLQQFRAMFMKRLMYNWRNWRGILVQILGLIISTFLLLKSHEFRYKKIRQMNLDEYGQTIVPFSIWGKSNLTSSLLTHLENMLKPGNHQLKEVQGDLLKYLEGNDECVHLCVIALSIKVVANRVNLTVLFNNEAYHSPSLSLTVLDNILFMSLSGSDASITVFNKPQPSPQRKEWPGSTDGKIVAFKIQLGMALLVSGFCILTVTERHNKTKHMQFLSGVSILVYWLSALVFDLIIFFISCCFLLVMFKYCKFDIYVTDYHILDTMLILTLFGWSAIPLTYLLSFLFSKSNSAYINLLVFCYLSGTLSLLMDTIIEARISTIMSNSTQTFLLNALLLFPMYNLGKCISEYTVIYRKKMLCIQQKNALKYLNCSNKHTKKNIYSLKKPMLGKYLIAMSIAGFVFLLLIFFWENISWKVKMFIHQHIYFGACKKYKPDIISKELSGTSEDNDVENERREILYQPEKFLNCPVLIKELTKIYFKSPLILAVKNISLAIQERACFGLLGFNGAGKTTTFQILTGENIPTAGDVFIDGISLTKNIVKVRSKIGYCPQFDALLEYMTGWEIMIMYARIWGISEHQIQPYVKKYLNSLDLESHANSLISTYSEGNKRRLSTAIATMGKPSVIFLDEPSTGMDPRARRLLWDTVIKIRESGKAIIITSHSMEECEALCTRLSIMVRGRLTCLGSPQYLKNKFGNIYILKAKVKSGETLDEFKNFITLTFPGSELQQENQGILNYCIPRKNNSWGKVFGILEKAKEQYNLEDYSISQITLDQVFLSFADQDRK"
        
        cutoff = 0.00001
        domains = pfam_tools.get_computed_pfam_domains(prot_seq, cutoff)
        
        self.assertEqual(2, len(domains))
        
        d = domains[0]
        self.assertEqual(('ABC_tran', 565, 685), (d.label, d.start, d.stop))
        
        d = domains[1]
        self.assertEqual(('ABC_tran', 1396, 1516), (d.label, d.start, d.stop))

    @patch('ptmworker.pfam_tools.get_computed_pfam_domains')
    def test_parse_or_query_domains_should_query_if_domains_empty(self, patch_query):
        d1 = pfam_tools.PFamDomain()
        d1.p_value = 0.01
        d1.start = 100
        d1.stop = 200
        d1.release = 23
        d1.label = 'Domain'

        prot = protein.Protein()
        patch_query.return_value = [d1]

        pfam_tools.parse_or_query_domains(prot, [])
        result = prot.domains

        self.assertEqual(1, len(result))
        self.assertEqual(0.01, result[0].p_value)
        self.assertEqual("COMPUTED PFAM", result[0].source)
        self.assertEqual(100, result[0].start)
        self.assertEqual(200, result[0].stop)
        self.assertEqual('Domain', result[0].label)
        self.assertEqual(23, result[0].version)


    def test_parse_or_query_domains_should_not_query_if_domains_provided(self):
        d1 = pfam_tools.PFamDomain()
        d1.p_value = 0.01
        d1.start = 100
        d1.stop = 200
        d1.release = 23
        d1.label = 'Domain'

        prot = protein.Protein()

        pfam_tools.parse_or_query_domains(prot, [d1])
        result = prot.domains

        self.assertEqual(1, len(result))
        self.assertEqual(0.01, result[0].p_value)
        self.assertEqual("PARSED PFAM", result[0].source)
        self.assertEqual(100, result[0].start)
        self.assertEqual(200, result[0].stop)
        self.assertEqual('Domain', result[0].label)
        self.assertEqual(23, result[0].version)

#    def test_get_computed_pfam_domains_2(self):
#        prot_seq = "MKKFFDSRREQGGSGLGSGSSGGGGSTSGLGSGYIGRVFGIGRQQVTVDEVLAEGGFAIVFLVRTSNGMKCALKRMFVNNEHDLQVCKREIQIMRDLSGHKNIVGYIDSSINNVSSGDVWEVLILMDFCRGGQVVNLMNQRLQTGFTENEVLQIFCDTCEAVARLHQCKTPIIHRDLKVENILLHDRGHYVLCDFGSATNKFQNPQTEGVNAVEDEIKKYTTLSYRAPEMVNLYSGKIITTKADIWALGCLLYKLCYFTLPFGESQVAICDGNFTIPDNSRYSQDMHCLIRYMLEPDPDKRPDIYQVSYFSFKLLKKECPIPNVQNSPIPAKLPEPVKASEAAAKKTQPKARLTDPIPTTETSIAPRQRPKAGQTQPNPGILPIQPALTPRKRATVQPPPQAAGSSNQPGLLASVPQPKPQAPPSQPLPQTQAKQPQAPPTPQQTPSTQAQGLPAQAQATPQHQQQLFLKQQQQQQQPPPAQQQPAGTFYQQQQAQTQQFQAVHPATQKPAIAQFPVVSQGGSQQQLMQNFYQQQQQQQQQQQQQQLATALHQQQLMTQQAALQQKPTMAAGQQPQPQPAAAPQPAPAQEPAIQAPVRQQPKVQTTPPPAVQGQKVGSLTPPSSPKTQRAGHRRILSDVTHSAVFGVPASKSTQLLQAAAAEASLNKSKSATTTPSGSPRTSQQNVYNPSEGSTWNPFDDDNFSKLTAEELLNKDFAKLGEGKHPEKLGGSAESLIPGFQSTQGDAFATTSFSAGTAEKRKGGQTVDSGLPLLSVSDPFIPLQVPDAPEKLIEGLKSPDTSLLLPDLLPMTDPFGSTSDAVIEKADVAVESLIPGLEPPVPQRLPSQTESVTSNRTDSLTGEDSLLDCSLLSNPTTDLLEEFAPTAISAPVHKAAEDSNLISGFDVPEGSDKVAEDEFDPIPVLITKNPQGGHSRNSSGSSESSLPNLARSLLLVDQLIDL"
#        
#        cutoff = 0.00001
#        domains = pfam_tools.get_computed_pfam_domains(prot_seq, cutoff)
#        
#        self.assertEqual(1, len(domains))
#        
#        testdom = domains[0]
#        self.assertEqual(46, testdom.start)
#        self.assertEqual(312, testdom.stop)
#        self.assertEqual(1, testdom.significant)
#        self.assertEqual("Pkinase", testdom.label)
#        self.assertEqual(4.5e-42, testdom.p_value)
#        self.assertEqual("26.0", testdom.release)
#        self.assertEqual("Domain", testdom.class_)
#        
#    
#    def test_get_computed_pfam_domains(self):
#        prot_seq = "QTTEIYFLSFGGCKAVTRAQHGQVLVRAPWLPSHRVLTWQRERDPKPSGVFFVCLLFFGTESHTAARLECSGPISAHCNLHLPGSSNSPAPASQLAGITGTHHHTRLIFVFLVERGFHHAGQAGLELLSPADPPASVSQSAGITGVSHCAWPELSPYRTSENPYKVRERPYKCIGCGKDLKGSSQLIQHQRIHSGKKLCECNECSKTCNQSSHFIRHHRIHSGEKPYECDKHGKAFRWSSDLVRHLKIHIGDKVLECSERGRIFNQSSGLMQHQRNHTES"
#        cutoff = 0.1
#        domains = pfam_tools.get_computed_pfam_domains(prot_seq, cutoff)
#        self.assertEqual(3, len(domains))
#        
#        testdom = domains[0]
#        self.assertEqual(185, testdom.start)
#        self.assertEqual(210, testdom.stop)
#        self.assertEqual(1, testdom.significant)
#        self.assertEqual("zf-H2C2_2", testdom.label)
#        self.assertEqual("PF13465.1", testdom.accession)
#        self.assertEqual(0.00093, testdom.p_value)
#        self.assertEqual("26.0", testdom.release)
#        self.assertEqual("Domain", testdom.class_)
        
