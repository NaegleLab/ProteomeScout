import os
from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.utils import downloadutils
from ptmscout.database import experiment, modifications
from ptmscout.config import settings

class TestUploadUtils(IntegrationTestCase):
    def test_export_experiment_28(self):
        exp = experiment.getExperimentById(28)

        new_ms_ids = {399217:'new_acc',
                      399218:'new_acc',
                      399219:'new_acc',
                      399220:'another'}

        headers, filename = downloadutils.experiment_to_tsv(exp, new_ms_ids)

        filepath = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, filename)

        exp_headers = ['accession', 'peptide', 'modification', 'run', \
                'data:time:0', 'data:time:5', 'data:time:10', 'data:time:30', \
                'stddev:time:0', 'stddev:time:10', 'stddev:time:30']
        self.assertEqual(exp_headers, headers)

        acc_map = {}
        run_names = set()
        i = 0
        with open(filepath, 'r') as tsv_file:
            for line in tsv_file:
                i+=1

                items = line.split('\t')

                acc = items[0]
                pep = items[1]

                if acc not in acc_map: 
                    acc_map[acc] = set()
                acc_map[acc].add( pep )

                run_names.add(items[3])

        for ms_id in new_ms_ids:
            self.assertIn(modifications.getMeasuredPeptideById(ms_id).peptide, acc_map[new_ms_ids[ms_id]])

        self.assertEqual(68*4+1, i)
        self.assertEqual(set(['run', '24H_EGF', '24H_HRG', 'P_EGF', 'P_HRG']), run_names)

        os.remove(filepath)

