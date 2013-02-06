import celery
import os
from tests.PTMScoutTestCase import IntegrationTestCase

@celery.task
def err_back(uuid, message, tmpfile):
    fn = open(tmpfile, 'w')

    result = celery.result.AsyncResult(uuid)
    exc = result.get(propagate=False)
    fn.write("%s\n%s\n" % ( str(exc), message ))

    fn.close()

@celery.task
def function(some_data):
    print "Trying but failing to do some stuff..."
    raise Exception(some_data)

class TestCase(IntegrationTestCase):
    def test_errbacks(self):
        tmpfile = 'errfile.tmp'
        res = function.apply_async(("some data",), link_error=err_back.s("an error occurred", tmpfile))
        r = res.get(propagate=False)

        self.assertTrue(os.path.exists('errfile.tmp'))

        fn = open(tmpfile, 'r')
        print fn.read()

        os.remove('errfile.tmp')
