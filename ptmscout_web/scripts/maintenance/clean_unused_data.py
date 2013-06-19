from scripts.DB_init import DatabaseInitialization
from ptmscout.config import settings
from ptmscout.database import DBSession, experiment, upload
import datetime
import time
import os
import re
import traceback
import shutil

DOWNLOAD_EXPIRATION_TIME = 60*60*24
SESSION_EXPIRATION_TIME = 60*60*24

def clean_old_files(now, path):
    for fn in os.listdir(path):
        fpath = os.path.join(path, fn)
        if os.path.isdir(fpath):
            clean_old_files( now, fpath )
            if os.listdir( fpath ) == []:
                print "Removing Dir:  " + fpath
                os.rmdir( fpath )
        else:
            ctime = os.path.getctime( fpath )
            if now - ctime > DOWNLOAD_EXPIRATION_TIME:
                print "Removing File: " + fpath
                os.remove( fpath )

if __name__ == "__main__":
    try:
        config_options = os.path.join(os.sep, 'data', 'ptmscout', 'ptmscout_web', 'production.ini')
        DatabaseInitialization.setUpClass(config_options)
        dbinit = DatabaseInitialization()
        dbinit.setUp()

        active_files = set()
        known_experiments = set()

        for exp in DBSession.query(experiment.Experiment):
            active_files.add(exp.dataset)
            known_experiments.add(exp.id)

        files_available = [ fn for fn in os.listdir(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path)) if fn.find('experiment')==0 ]

        active_sessions = []
        inactive_sessions = []
        now = datetime.datetime.now()

        sessions_by_exp = {}
        for session in DBSession.query(upload.Session):
            delta = now - session.date
            if session.experiment_id in known_experiments:
                value = sessions_by_exp.get(session.experiment_id, [])
                value.append((delta, session))
                sessions_by_exp[session.experiment_id] = value
            elif delta.total_seconds() >= SESSION_EXPIRATION_TIME:
                inactive_sessions.append(session)
            else:
                active_files.add(session.data_file)
                active_sessions.append(session)

        for exp_id in sessions_by_exp:
            sorted_sessions = sorted(sessions_by_exp[exp_id], key=lambda item: item[0])
            active_sessions.append(sorted_sessions[0][1])
            active_files.add(sorted_sessions[0][1].data_file)
            for delta, session in sorted_sessions[1:]:
                inactive_sessions.append(session)


        print "Keeping %d active or referenced files" % (len(active_files))
        print "Keeping %d active or referenced sessions" % (len(active_sessions))
        print "Deleting %d outdated sessions" % (len(inactive_sessions))

        print "Active/Referenced Sessions:"
        for session in active_sessions:
            print session.id, session.experiment_id, session.date

        print "Deleted Sessions:"
        for session in inactive_sessions:
            print session.id, session.experiment_id, session.date
            session.delete()

        print "Deleted inactive files:"
        for fn in files_available:
            if fn not in active_files:
                print fn
                os.remove(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, fn))

        print "Deleted inactive experiment directories:"
        for dn in os.listdir(os.path.join(settings.ptmscout_path, settings.experiment_data_file_path)):
            dp = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, dn)
            m = re.match(r'^e([0-9]+)$', dn)
            if os.path.isdir(dp) and m:
                exp_id = int(m.group(1))
                if exp_id not in known_experiments:
                    print dn
                    shutil.rmtree(dp)

        now = time.time()
        print "Cleaning exported MCAM and experiment files"
        clean_old_files( now, os.path.join(settings.ptmscout_path, settings.annotation_export_file_path) )
        clean_old_files( now, os.path.join(settings.ptmscout_path, settings.mcam_file_path) )

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()
