from scripts.DB_init import DatabaseInitialization
from ptmscout.config import settings
from ptmscout.database import DBSession, experiment, upload
from paste.deploy.loadwsgi import appconfig
import datetime
import sys, os
import traceback

SESSION_EXPIRATION_TIME = 86400

if __name__ == "__main__":
    try:
        config_options = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'production.ini'))
        
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

        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()