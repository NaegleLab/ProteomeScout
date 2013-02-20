import os
from tests.PTMScoutTestCase import IntegrationTestCase
from ptmscout.utils import downloadutils
from ptmscout.database import experiment
from ptmscout.config import settings

class TestUploadUtils(IntegrationTestCase):
    def test_export_experiment_26(self):
        exp = experiment.getExperimentById(28)
        headers, filename = downloadutils.experiment_to_tsv(exp)

        filepath = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, filename)

        exp_headers = ['accession', 'peptide', 'modification', 'run', \
                'data:time:0', 'data:time:5', 'data:time:10', 'data:time:30', \
                'stddev:time:0', 'stddev:time:10', 'stddev:time:30']
        self.assertEqual(exp_headers, headers)

        run_names = set()
        i = 0
        with open(filepath, 'r') as tsv_file:
            for line in tsv_file:
                i+=1

                items = line.split('\t')
                run_names.add(items[3])

        self.assertEqual(68*4+1, i)
        self.assertEqual(set(['run', '24H_EGF', '24H_HRG', 'P_EGF', 'P_HRG']), run_names)

        os.remove(filepath)

