from scripts.DB_init import DatabaseInitialization
from ptmscout.database import jobs
import traceback
import os, sys

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

        running, error = [], []
        for j in dbinit.session.query(jobs.Job):
            if j.status != 'finished' and j.status !='error':
                running.append(j)
            if j.status == 'error':
                error.append(j)


        print "\n-----------\nRunning jobs:\n-----------"
        print_jobs(running)

        print "\n-----------\nError jobs:\n-----------"
        print_jobs(error)
    finally:
        dbinit.close()
