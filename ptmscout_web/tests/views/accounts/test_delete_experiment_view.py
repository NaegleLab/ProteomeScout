from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.config import strings
from ptmscout.database import experiment

class IntegrationTestDeleteExperimentView(IntegrationTestCase):
    def test_view_integration(self):
        exp = experiment.getExperimentById(26)
        exp.type='dataset'
        exp.saveExperiment()
        
        self.bot.acquire_experiments([26])
        
        result = self.ptmscoutapp.get("/experiments/26/delete", status=200)
        
        result.mustcontain(strings.delete_experiment_confirm_message)
        result = result.forms[0].submit()
        
        result.mustcontain(strings.delete_experiment_success_message)
        
        try:
            experiment.getExperimentById(26, self.bot.user)
        except experiment.NoSuchExperiment:
            pass
        else:
            self.fail("Expected no such experiment exception")