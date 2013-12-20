from scripts.DB_init import DatabaseInitialization
from ptmscout.database import jobs
from ptmscout.utils import mail
from ptmscout.config import strings, settings
import traceback
import os, sys

def notify_job_failed(job, exc, stack_trace):
    job.fail(stack_trace)
    job.save()
    
    subject = strings.job_failed_subject
    message = strings.job_failed_message % (job.name, job.stage, "Exception: " + str(exc))
    
    mail.celery_send_mail([job.user.email, settings.adminEmail], subject, message)
 

def print_jobs(job_list):
    for j in job_list:
        tp = j.type
        name = j.name
        if len(name) > 30:
            name = name[:27] + "..."
        uname = j.user.name
        created = j.created
        print "%-20s %-20s %-30s %20s" % (str(created), uname, name, tp)

if __name__ == "__main__":
    try:
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        running  = []
        for j in dbinit.session.query(jobs.Job):
            if j.status != 'finished' and j.status !='error':
                running.append(j)

        if len(running)==0:
            print "No running jobs to be canceled!"
        else:
            print "The following jobs will be cancelled:"
            print "\n-----------\nRunning jobs:\n-----------"
            print_jobs(running)

            print "\n\n"
            val = ""
            while val not in ['yes','no','n','y']:
                val = raw_input("Continue? (yes/no) ").lower()

            if val in ['yes','y']:
                for j in running:
                    notify_job_failed(j, "Server is shutting down, please try again later.", "Server shutdown!")

            dbinit.session.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
