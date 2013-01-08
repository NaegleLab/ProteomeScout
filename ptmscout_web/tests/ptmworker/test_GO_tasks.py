from ptmworker import GO_tasks
from tests.PTMScoutTestCase import IntegrationTestCase
from tests.views.mocking import createMockExperiment, createMockError
import pickle
from mock import patch
from ptmscout.config import settings, strings
from ptmscout.database import experiment, protein

class GOTasksTest(IntegrationTestCase):
    @patch('transaction.commit')
    def test_import_go_terms(self, patch_commit):
        def save_db():
            from ptmscout.database import DBSession
            DBSession.flush()

        def show_go_terms(prots):
            for pid in prots:
                p = protein.getProteinById(pid)
                print [ term.GO_term.GO for term in p.GO_terms ]

        protein_map = {"some protein results"}
        new_protein_ids = {"Q9NQG7": 9, "Q9WTS6": 11, "P49848":24, "O88384":105}

        patch_commit.side_effect = save_db
        protein_result = protein_map, new_protein_ids

        result = GO_tasks.import_go_terms.apply_async((protein_result, 26))

        value = result.get(propagate=False)
        if isinstance(value, Exception):
            print result.traceback
            self.fail()

        self.assertEqual(protein_map, value)

        p = protein.getProteinById(11)
        self.assertEqual(['GO:0005887', 'GO:0016020', 'GO:0046982',
                'GO:0007156', 'GO:0010976', 'GO:0042803'], [ term.GO_term.GO
                    for term in p.GO_terms])

        p = protein.getProteinById(9)
        self.assertIn('GO:0016020', [ term.GO_term.GO for term in p.GO_terms ])


        exp = experiment.getExperimentById(26, None, False, False)
        self.assertEqual('GO terms', exp.loading_stage)

        self.assertTrue(patch_commit.called)
