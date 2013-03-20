from tests.DBTestCase import DBTestCase
import ptmscout.database.experiment as dbexperiment
from ptmscout.database.experiment import NoSuchExperiment, Experiment,\
    getExperimentTree, ExperimentAccessForbidden, ExperimentNotAvailable,\
    searchExperiments
from mock import patch
from ptmscout.database import user, permissions

class ExperimentTestCase(DBTestCase):

    def test_experiment_searchExperiments_by_condition(self):
        cnt, experiments = searchExperiments(text_search='HER2', conditions={'drug':['dasatinib']})
        
        self.assertEqual(1, cnt)
        self.assertEqual(1, len(experiments))
        self.assertEqual('Effects of HER2 overexpression on cell signaling networks governing proliferation and migration.', experiments[0].name)

    def test_experiment_searchExperiments_by_name(self):
        cnt, experiments = searchExperiments(text_search='EGFR')

        self.assertEqual(1, cnt)
        self.assertEqual(1, len(experiments))
        self.assertEqual('Quantitative analysis of EGFRvIII cellular signaling networks reveals a combinatorial therapeutic strategy for glioblastoma.', experiments[0].name)

    def format_multiline(self, string):
        return string.replace("    ","").replace("\n"," ")

    def test_experiment_getCitationString_should_handle_multiple_page_range_types(self):
        exp = dbexperiment.getExperimentById(26, None)
        exp.page_start = "571"

        exp.page_end = "571"
        exp_cite = self.format_multiline("""Mol Cell Proteomics. 2005-september. Vol 4. 571.""")

        self.assertEqual(exp_cite, exp.getCitationString())

        exp.page_end = None
        exp_cite = self.format_multiline("""Mol Cell Proteomics. 2005-september. Vol 4. 571.""")

        self.assertEqual(exp_cite, exp.getCitationString())

        exp.page_end = "5712"
        exp_cite = self.format_multiline("""Mol Cell Proteomics. 2005-september. Vol 4. 571-5712.""")

        self.assertEqual(exp_cite, exp.getCitationString())

    def test_experiment_getLongCitationString_should_handle_multiple_page_range_types(self):
        exp = dbexperiment.getExperimentById(26, None)
        exp.page_start = "571"

        exp.page_end = "571"
        exp_cite = self.format_multiline("""Yi Zhang, Alejandro Wolf-Yadlin, Phillip L Ross, Darryl J Pappin, John Rush,
        Douglas A Lauffenburger, Forest M White. <b>Mol Cell Proteomics</b>. 2005-september. Vol 4. 571.""")

        self.assertEqual(exp_cite, exp.getLongCitationString())

        exp.page_end = None
        exp_cite = self.format_multiline("""Yi Zhang, Alejandro Wolf-Yadlin, Phillip L Ross, Darryl J Pappin, John Rush,
        Douglas A Lauffenburger, Forest M White. <b>Mol Cell Proteomics</b>. 2005-september. Vol 4. 571.""")

        self.assertEqual(exp_cite, exp.getLongCitationString())

        exp.page_end = "5712"
        exp_cite = self.format_multiline("""Yi Zhang, Alejandro Wolf-Yadlin, Phillip L Ross, Darryl J Pappin, John Rush,
        Douglas A Lauffenburger, Forest M White. <b>Mol Cell Proteomics</b>. 2005-september. Vol 4. 571-5712.""")

        self.assertEqual(exp_cite, exp.getLongCitationString())

    def test_experiment_copy_should_copy_all_data(self):
        exp = dbexperiment.getExperimentById(26, None)
        
        nexp = dbexperiment.Experiment()
        nexp.copyData(exp)
        
        self.assertEqual(exp.name, nexp.name)
        self.assertEqual(exp.author, nexp.author)
        self.assertEqual(exp.date, nexp.date)

        self.assertEqual(exp.description, nexp.description)
        self.assertEqual(exp.contact, nexp.contact)

        self.assertEqual(exp.PMID, nexp.PMID)
        self.assertEqual(exp.URL, nexp.URL)
        self.assertEqual(exp.published, nexp.published)

        self.assertEqual(exp.ambiguity, nexp.ambiguity)
        self.assertEqual(exp.export, nexp.export)
        self.assertEqual(exp.experiment_id, nexp.experiment_id)
        
        self.assertEqual(exp.dataset, nexp.dataset)
        
        self.assertEqual(exp.volume, nexp.volume)
        self.assertEqual(exp.page_start, nexp.page_start)
        self.assertEqual(exp.page_end, nexp.page_end)
        self.assertEqual(exp.journal, nexp.journal)
        self.assertEqual(exp.publication_year, nexp.publication_year)
        self.assertEqual(exp.publication_month, nexp.publication_month)

        self.assertEqual(exp.public, nexp.public)

        self.assertEqual(exp.submitter_id, nexp.submitter_id)

        self.assertEqual([], nexp.errors)
        self.assertEqual([], nexp.conditions)
        
        
        
    
    def test_getExperimentById_should_succeed_on_existing_experiment(self):
        experiment = dbexperiment.getExperimentById(1, None)
        
        self.assertEqual("PhosphoSite: A bioinformatics resource dedicated to physiological protein phosphorylation.", experiment.name)
        self.assertEqual("http://www.phosphosite.org/", experiment.URL)
    
    def test_getExperimentById_should_throw_exception_when_experiment_is_currently_uploading(self):
        try:
            exp = dbexperiment.getExperimentById(1, None)
            exp.public = 0
            exp.job.status = 'loading'
            exp.saveExperiment()
            
            ptmuser = user.getUserById(1)
            
            perm = permissions.Permission(exp)
            perm.experiment_id = exp.id
            perm.user_id = ptmuser.id
            
            ptmuser.permissions.append(perm)
            
            dbexperiment.getExperimentById(1, ptmuser)
            
        except ExperimentNotAvailable, f:
            self.assertEqual(1, f.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception ExperimentNotAvailable was not thrown")
    
    def test_getExperimentById_should_throw_exception_when_experiment_is_private(self):
        try:
            exp = dbexperiment.getExperimentById(1, None)
            exp.public = 0
            exp.saveExperiment()
            dbexperiment.getExperimentById(1, None)
            
        except ExperimentAccessForbidden, f:
            self.assertEqual(1, f.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchExperiment was not thrown")
        
    
    def test_getExperimentById_should_throw_exception_when_no_such_experiment(self):
        try:
            dbexperiment.getExperimentById(100000, None)
        except NoSuchExperiment, nse:
            self.assertEqual(100000, nse.eid)
        except Exception, e:
            self.fail("Unexpected exception: " + str(e))
        else:
            self.fail("Expected exception NoSuchExperiment was not thrown")
            
    def test_getAllExperiments_should_return_all_experiments_in_database_to_which_access_is_allowed(self):
        exp = dbexperiment.getExperimentById(1, None)
        exp.public = 0
        exp.saveExperiment()
        
        exp2 = dbexperiment.getExperimentById(26, None)
        exp2.job.status = 'configuration'
        exp2.saveExperiment()
        
        experiments = dbexperiment.getAllExperiments(None)
        
        self.assertTrue(len(experiments) > 0)
        
        self.assertFalse( exp.id in [ e.id for e in experiments ] )
        self.assertFalse( exp2.id in [ e.id for e in experiments ] )
        
    @patch('ptmscout.database.experiment.getAllExperiments')
    def test_getExperimentTree_should_build_experiment_tree(self, patch_getAllExperiments):
        experiments = [Experiment(), Experiment(), Experiment(), Experiment(), Experiment(), Experiment()]
        
        experiments[0].id = 1
        experiments[1].id = 2
        experiments[2].id = 3
        experiments[3].id = 4
        experiments[4].id = 5
        experiments[5].id = 6
        
        experiments[0].experiment_id = None
        experiments[1].experiment_id = None
        experiments[2].experiment_id = 1
        experiments[3].experiment_id = 1
        experiments[4].experiment_id = None
        experiments[5].experiment_id = 2
        
        patch_getAllExperiments.return_value = experiments
        
        current_user = 0
        experiment_tree = getExperimentTree(current_user)
        
        patch_getAllExperiments.assert_called_with(current_user)
        self.assertEqual([experiments[0], experiments[1], experiments[4]], experiment_tree)
        self.assertEqual([experiments[2], experiments[3]], experiments[0].children)
        self.assertEqual([experiments[5]], experiments[1].children)
        
