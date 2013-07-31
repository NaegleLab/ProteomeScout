from mock import patch, Mock
from ptmscout.config import strings, settings
from pyramid.testing import DummyRequest
from tests.PTMScoutTestCase import UnitTestCase, IntegrationTestCase
from ptmscout.views.protein.structure_view import protein_structure_viewer, \
            format_protein_modifications, format_protein_domains
from tests.views.mocking import createMockProtein, createMockUser, \
    createMockMeasurement, createMockExperiment, createMockPeptide,\
    createMockPTM, createMockPeptideModification, createMockDomain,\
    createMockData, createMockRegion
import json, base64

class TestProteinStructureViewIntegration(IntegrationTestCase):
    def test_default_page_is_structure_view(self):
        result = self.ptmscoutapp.get('/proteins/35546')
        d = result.pyquery
        encoded_data = d(".protein_viewer .data").text()

        protein_data = json.loads( base64.b64decode( encoded_data ))

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




    def test_integration(self):
        result = self.ptmscoutapp.get('/proteins/35546/sequence')
        d = result.pyquery
        encoded_data = d(".protein_viewer .data").text()

        protein_data = json.loads( base64.b64decode( encoded_data ))

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
        request = DummyRequest()
        request.route_url = Mock()
        request.route_url.return_value = 'some_experiment_page'
        prot = createMockProtein()

        exp1 = createMockExperiment(13)
        exp2 = createMockExperiment(16)

        exp1.type = 'experiment'
        exp2.type = 'compendia'
        exp2.URL = 'http://compendia.com'

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

        createMockPeptideModification(ms1, pep1, ptm1)
        createMockPeptideModification(ms1, pep2, ptm1)
        createMockPeptideModification(ms2, pep2, ptm2)
        createMockPeptideModification(ms2, pep3, ptm1)

        region1 = createMockRegion(pid=prot.id, label='kinase', start=190, stop=220)
        prot.regions.append(region1)

        experiments, mod_types, result = format_protein_modifications(request, prot, measured_peps)

        self.maxDiff=None

        exp_result = {  200: {'mods':{}, 'residue':pep1.site_type, 'domain':None, 'peptide': pep1.pep_aligned, 'regions':['kinase']},
                        210: {'mods':{}, 'residue':pep2.site_type, 'domain':'d1', 'peptide': pep2.pep_aligned, 'regions':['kinase']},
                        300: {'mods':{}, 'residue':pep3.site_type, 'domain':None, 'peptide': pep3.pep_aligned, 'regions':[] }}

        exp_result[200]['mods'][ptm1.name] = [ {'MS': ms1.id, 'experiment_url': 'some_experiment_page', 'experiment':13, 'has_data': True} ]
        exp_result[210]['mods'][ptm1.name] = [ {'MS': ms1.id, 'experiment_url': 'some_experiment_page', 'experiment':13, 'has_data': True} ]
        exp_result[210]['mods'][ptm2.name] = [ {'MS': ms2.id, 'experiment_url': 'http://compendia.com', 'experiment':16, 'has_data': False} ]
        exp_result[300]['mods'][ptm1.name] = [ {'MS': ms2.id, 'experiment_url': 'http://compendia.com', 'experiment':16, 'has_data': False} ]

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

    @patch('ptmscout.views.protein.structure_view.format_scansite_predictions')
    @patch('ptmscout.views.protein.structure_view.format_protein_mutations')
    @patch('ptmscout.views.protein.structure_view.format_protein_regions')
    @patch('ptmscout.views.protein.structure_view.format_protein_domains')
    @patch('ptmscout.views.protein.structure_view.format_protein_modifications')
    @patch('ptmscout.database.modifications.getMeasuredPeptidesByProtein')
    @patch('ptmscout.database.protein.getProteinById')
    def test_protein_structure_view(self, patch_getProtein, patch_getMods,
            patch_formatMods, patch_formatDomains, patch_formatRegions,
            patch_formatMutations, patch_formatScansite):

        request = DummyRequest()
        request.user = createMockUser()
        request.matchdict['id'] = '35546'
        request.GET['experiment_id']='1302'

        prot = createMockProtein()
        prot.id = 35546
        patch_getProtein.return_value = prot

        formatted_domains = "some formatted domains"
        formatted_exps = "some experiments"
        formatted_mod_types = "some mod types"
        formatted_mods = "some formatted modifications"
        formatted_muts = "some formatted mutations"
        formatted_regions = "some formatted regions"
        formatted_scansite = "some formatted scansite data"

        exp_mods = ["some", "mods"]
        patch_getMods.return_value = exp_mods
        patch_formatDomains.return_value = formatted_domains
        patch_formatMods.return_value = (formatted_exps, formatted_mod_types, formatted_mods)
        patch_formatRegions.return_value = formatted_regions
        patch_formatMutations.return_value = formatted_muts
        patch_formatScansite.return_value = formatted_scansite

        result = protein_structure_viewer(None, request)

        patch_getMods.assert_called_once_with(prot.id, request.user)
        patch_getProtein.assert_called_once_with(35546)
        patch_formatMods.assert_called_once_with(request, prot, exp_mods)
        patch_formatDomains.assert_called_once_with(prot)
        patch_formatRegions.assert_called_once_with(prot)
        patch_formatMutations.assert_called_once_with(prot)
        patch_formatScansite.assert_called_once_with(prot)

        self.assertEqual(strings.protein_structure_page_title, result['pageTitle'])
        self.assertEqual(prot, result['protein'])
        decoded_data = json.loads(base64.b64decode(result['data']))

        exp_data = {'seq': prot.sequence,
                    'domains': formatted_domains,
                    'mods': formatted_mods,
                    'exps': formatted_exps,
                    'mod_types': formatted_mod_types,
                    'mutations': formatted_muts,
                    'scansite': formatted_scansite,
                    'regions': formatted_regions,
                    'pfam_url': settings.pfam_family_url,
                    'protein_data_url': "%s/proteins/%s/data" % (request.application_url, prot.id),
                    'experiment':1302
                }

        exp_tracks = [
                    'PFam Domains',
                    'PTMs',
                    "Activation Loops",
                    "Uniprot Domains",
                    "Uniprot Structure",
                    "Uniprot Binding Sites",
                    "Uniprot Macrostructure",
                    "Uniprot Topology",
                    "Entrez Domains",
                    'Mutations',
                    'Scansite',
                    ]

        self.assertEqual(exp_data, decoded_data)
        self.assertEqual(formatted_exps, result['experiments'])
        self.assertEqual(formatted_mod_types, result['mod_types'])
        self.assertEqual(exp_tracks, result['tracks'])
        self.assertEqual("/".join([settings.documentationUrl, settings.proteinViewerHelp]), result['protein_viewer_help_page'])
