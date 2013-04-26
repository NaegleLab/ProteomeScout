from tests.PTMScoutTestCase import IntegrationTestCase, UnitTestCase
from ptmscout.database import upload
from ptmscout.config import strings

class IntegrationTestConfigureAnnotationUploadView(IntegrationTestCase):
    def test_view_should_show_errors_when_missing_columns(self):
        session = upload.Session()
        session.user_id = self.bot.user.id
        session.change_name=''
        session.change_description = ""
        session.data_file = 'test/experiment.28.test.annotations.tsv'
        session.resource_type = 'annotations'
        session.load_type = 'new'
        session.stage='config'
        session.save()
        
        result = self.ptmscoutapp.get("/experiments/28/annotate/%d/configure" % (session.id), status=200)
        result.form.set('column_0_type', 'none')
        result.form.set('column_34_type', 'none')
        result.form.set('column_35_type', 'none')
        result.form.set('column_36_type', 'none')
        result = result.form.submit()
        
        result.mustcontain(strings.experiment_upload_error_no_column_assignment % ('MS_id'))
        result.mustcontain(strings.experiment_upload_error_no_annotations)
        
    def test_view_should_fail_on_bad_value_redirect_on_override(self):
        session = upload.Session()
        session.user_id = self.bot.user.id
        session.change_name=''
        session.change_description = ""
        session.data_file = 'test/experiment.28.test.annotations.tsv'
        session.resource_type = 'annotations'
        session.load_type = 'new'
        session.stage='config'
        session.save()
        
        # get the form
        result = self.ptmscoutapp.get("/experiments/28/annotate/%d/configure" % (session.id), status=200)
        
        # submit with bad value and no override should fail
        result = result.form.submit()
        result.mustcontain(strings.experiment_upload_warning_columns_values_should_be % ('integer'))
        result.form.set('override', False)
        result = result.form.submit(status=200)
        
        # submit with bad value + override should work
        result.mustcontain(strings.experiment_upload_warning_columns_values_should_be % ('integer'))
        result.form.set('override', True)
        result = result.form.submit(status=302)
        
        self.assertEqual('MS_id', session.columns[0].type)
        self.assertEqual('cluster', session.columns[34].type)
        self.assertEqual('numeric', session.columns[35].type)
        self.assertEqual('nominative', session.columns[36].type)
