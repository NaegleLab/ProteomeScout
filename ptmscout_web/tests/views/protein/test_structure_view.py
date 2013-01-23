from mock import patch
from ptmscout.config import strings, settings
from pyramid.testing import DummyRequest
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
from ptmscout.views.protein.structure_view import protein_structure_viewer, \
            format_protein_modifications, format_protein_domains
from tests.views.mocking import createMockProtein, createMockUser, \
    createMockMeasurement, createMockExperiment, createMockPeptide,\
    createMockPTM, createMockPeptideModification, createMockDomain,\
    createMockData
import json, base64

class TestProteinStructureViewIntegration(IntegrationTestCase):
    def test_integration(self):
        result = self.ptmscoutapp.get('/proteins/35546/structure')
        d = result.pyquery
        protein_data = json.loads( base64.b64decode( d(".protein_viewer .data").text() ))
        self.maxDiff=None

        self.assertEqual(
            {'25': "Quantitative analysis of EGFRvIII cellular signaling networks reveals a combinatorial therapeutic strategy for glioblastoma.",
             '26': "Time-resolved mass spectrometry of tyrosine phosphorylation sites in the epidermal growth factor receptor signaling network reveals dynamic modules.",
             '28': "Effects of HER2 overexpression on cell signaling networks governing proliferation and migration." },
            protein_data['exps'])

        self.assertEqual(2, len(protein_data['domains']))
        self.assertEqual(['Phosphotyrosine'], protein_data['mod_types'])
        self.assertEqual(3, len(protein_data['mods']))
        self.assertEqual(1, len( protein_data['mods']['518']['mods']['Phosphotyrosine']))
        self.assertEqual(3, len( protein_data['mods']['857']['mods']['Phosphotyrosine']))
        self.assertEqual(1, len( protein_data['mods']['858']['mods']['Phosphotyrosine']))
        self.assertEqual(protein_data['seq'], "MQPEEGTGWLLELLSEVQLQQYFLRLRDDLNVTRLSHFEYVKNEDLEKIGMGRPGQRRLWEAVKRRKALCKRKSWMSKVFSGKRLEAEFPPHHSQSTFRKTSPAPGGPAGEGPLQSLTCLIGEKDLRLLEKLGDGSFGVVRRGEWDAPSGKTVSVAVKCLKPDVLSQPEAMDDFIREVNAMHSLDHRNLIRLYGVVLTPPMKMVTELAPLGSLLDRLRKHQGHFLLGTLSRYAVQVAEGMGYLESKRFIHRDLAARNLLLATRDLVKIGDFGLMRALPQNDDHYVMQEHRKVPFAWCAPESLKTRTFSHASDTWMFGVTLWEMFTYGQEPWIGLNGSQILHKIDKEGERLPRPEDCPQDIYNVMVQCWAHKPEDRPTFVALRDFLLEAQPTDMRALQDFEEPDKLHIQMNDVITVIEGRAENYWWRGQNTRTLCVGPFPRNVVTSVAGLSAQDISQPLQNSFIHTGHGDSDPRHCWGFPDRIDELYLGNPMDPPDLLSVELSTSRPPQHLGGVKKPTYDPVSEDQDPLSSDFKRLGLRKPGLPRGLWLAKPSARVPGTKASRGSGAEVTLIDFGEEPVVPALRPCPPSLAQLAMDACSLLDETPPQSPTRALPRPLHPTPVVDWDARPLPPPPAYDDVAQDEDDFEICSINSTLVGAGVPAGPSQGQTNYAFVPEQARPPPPLEDNLFLPPQGGGKPPSSAQTAEIFQALQQECMRQLQAPGSPAPSPSPGGDDKPQVPPRVPIPPRPTRPHVQLSPAPPGEEETSQWPGPASPPRVPPREPLSPQGSRTPSPLVPPGSSPLPPRLSSSPGKTMPTTQSFASDPKYATPQVIQAPGAGGPCILPIVRDGKKVSSTHYYLLPERPSYLERYQRFLREAQSPEEPTPLPVPLLLPPPSTPAPAAPTATVRPMPQAALDPKANFSTNNSNPGARPPPPRATARLPQRGCPGDGPEAGRPADKIQMAMVHGVTTEECQAALQCHGWSVQRAAQYLKVEQLFGLGLRPRGECHKVLEMFDWNLEQAGCHLLGSWGPAHHKR")


class TestProteinStructureViews(UnitTestCase):

    def test_format_protein_modifications(self):
        prot = createMockProtein()

        exp1 = createMockExperiment(13)
        exp2 = createMockExperiment(16)

        exp1.export = 1
        exp2.export = 0

        ms1 = createMockMeasurement(prot.id, 13)
        ms1.experiment = exp1
        ms2 = createMockMeasurement(prot.id, 16)
        ms2.experiment = exp2
        measured_peps = [ms1, ms2]

        ms1.data = createMockData(6, 'avg', ms1.id)

        ptm1 = createMockPTM()
        ptm2 = createMockPTM()

        pep1 = createMockPeptide(prot.id, 200)
        pep2 = createMockPeptide(prot.id, 210)
        pep3 = createMockPeptide(prot.id, 300)

        pep1.protein_domain = None
        pep2.protein_domain = createMockDomain(prot.id, 'd1')
        pep3.protein_domain = None

        pepmod1 = createMockPeptideModification(ms1, pep1, ptm1)
        pepmod2 = createMockPeptideModification(ms1, pep2, ptm1)
        pepmod3 = createMockPeptideModification(ms2, pep2, ptm2)
        pepmod4 = createMockPeptideModification(ms2, pep3, ptm1)

        experiments, mod_types, result = format_protein_modifications(prot, measured_peps)

        exp_result = {  200: {'mods':{}, 'residue':pep1.site_type, 'domain':None, 'peptide': pep1.pep_aligned },
                        210: {'mods':{}, 'residue':pep2.site_type, 'domain':'d1', 'peptide': pep2.pep_aligned },
                        300: {'mods':{}, 'residue':pep3.site_type, 'domain':None, 'peptide': pep3.pep_aligned }}

        exp_result[200]['mods'][ptm1.name] = [ {'MS': ms1.id, 'experiment':13, 'has_data': True, 'exported':True} ]
        exp_result[210]['mods'][ptm1.name] = [ {'MS': ms1.id, 'experiment':13, 'has_data': True, 'exported':True} ]
        exp_result[210]['mods'][ptm2.name] = [ {'MS': ms2.id, 'experiment':16, 'has_data': False, 'exported':False} ]
        exp_result[300]['mods'][ptm1.name] = [ {'MS': ms2.id, 'experiment':16, 'has_data': False, 'exported':False} ]

        self.assertEqual(exp_result, result)
        self.assertEqual(sorted([ptm1.name, ptm2.name]), mod_types)
        self.assertEqual({13: exp1.name, 16: exp2.name}, experiments)

    def test_format_protein_domains(self):
        prot = createMockProtein()
        d1 = createMockDomain(prot.id, 'd1')
        d2 = createMockDomain(prot.id, 'd2')

        prot.domains = [d1,d2]

        formatted_domains = format_protein_domains(prot)

        exp_domains = \
        [{'label': d1.label, 'start': d1.start, 'stop': d1.stop, 'source': d1.source}, \
         {'label': d2.label, 'start': d2.start, 'stop': d2.stop, 'source': d2.source}]

        self.assertEqual(sorted(exp_domains, key=lambda d: d['start']), formatted_domains)

    @patch('ptmscout.views.protein.structure_view.format_protein_domains')
    @patch('ptmscout.views.protein.structure_view.format_protein_modifications')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_structure_view(self, patch_getProtein, patch_getMods,
            patch_formatMods, patch_formatDomains):
        request = DummyRequest()
        request.user = createMockUser()
        request.matchdict['id'] = '35546'
        request.GET['experiment_id']=1302

        prot = createMockProtein()
        prot.id = 35546
        patch_getProtein.return_value = prot

        formatted_domains = "some formatted domains"
        formatted_exps = "some experiments"
        formatted_mod_types = "some mod types"
        formatted_mods = "some formatted modifications"
        exp_mods = ["some", "mods"]
        patch_getMods.return_value = exp_mods
        patch_formatDomains.return_value = formatted_domains
        patch_formatMods.return_value = (formatted_exps, formatted_mod_types, formatted_mods)

        result = protein_structure_viewer(request)

        patch_getMods.assert_called_once_with(prot.id, request.user)
        patch_getProtein.assert_called_once_with(35546)
        patch_formatMods.assert_called_once_with(prot, exp_mods)
        patch_formatDomains.assert_called_once_with(prot)

        self.assertEqual(strings.protein_structure_page_title, result['pageTitle'])
        self.assertEqual(prot, result['protein'])
        decoded_data = json.loads(base64.b64decode(result['data']))

        exp_data = {'seq': prot.sequence,
                    'domains': formatted_domains,
                    'mods': formatted_mods,
                    'exps': formatted_exps,
                    'mod_types': formatted_mod_types,
                    'pfam_url': settings.pfam_family_url,
                    'protein_data_url': "%s/proteins/%s/data" % (request.application_url, prot.id),
                    'experiment_url': "%s/experiments" % (request.application_url),
                    'experiment':1302
                }

        self.assertEqual(exp_data, decoded_data)
        self.assertEqual(formatted_exps, result['experiments'])
        self.assertEqual(formatted_mod_types, result['mod_types'])
        self.assertEqual(['Domains', 'PTMs'], result['tracks'])
