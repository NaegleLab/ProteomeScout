import time
import os

def get_session(match_field, resource_type):
    def do_wrap(fn):
        from ptmscout.database import upload
        
        def wrapper(*args):
            request = args[1]
            sid = int(request.matchdict[match_field])
            session = upload.getSessionById(sid, user=request.user)
            
            if session.resource_type != resource_type:
                raise upload.NoSuchSession(sid) 
            
            nargs = tuple( list(args) + [session] )
            
            return fn(*nargs)
        return wrapper
    return do_wrap


def get_experiment(match_field, types=set(['experiment','dataset','compendia'])):
    def do_wrap(fn):
        from ptmscout.database import experiment
        
        def wrapper(*args):
            request = args[1]
            eid = int(request.matchdict[match_field])
            exp = experiment.getExperimentById(eid, user=request.user)
            
            if exp.type not in types:
                raise experiment.NoSuchExperiment(eid) 
            
            nargs = tuple( list(args) + [exp] )
            
            return fn(*nargs)
        
        return wrapper
    return do_wrap
            

def page_query(start=0, yield_per=1000):
    """
        Decorator to page queries
    """
    def wrapper(query_generator):
        offset = start
        while True:
            r = False
            for elem in query_generator(yield_per, offset):
                r = True
                yield elem
            offset += yield_per

            if not r: break
    return wrapper

def profile(fn):
    def do_profile(*args, **kwargs):
        t = time.clock()
        result = fn(*args, **kwargs)

        elapsed = time.clock() - t
        fn.profile_calls.append(elapsed)
        
        return result
    
    fn.profile_calls = []
    
    return do_profile
    

def pushdir(new_dir):
    def wrap(fn):
        def push_stack(*args, **kwargs):
            cwd = os.getcwd()
            os.chdir(new_dir)
            result = fn(*args, **kwargs)
            os.chdir(cwd)
            
            return result
        
        push_stack.__name__ = fn.__name__
        return push_stack
    
    return wrap

def rate_limit(rate=None):
    def wrap(fn):
        def rate_limited_task(*args, **kwargs):
            now = time.clock()
            diff = now - rate_limited_task.last_call_time

            if diff < rate_limited_task.seconds_per_task:
                time.sleep(rate_limited_task.seconds_per_task - diff)

            rval = fn(*args, **kwargs)

            rate_limited_task.last_call_time = now
            return rval

        if rate != None:
            rate_limited_task.seconds_per_task = 1.0 / float(rate)
        else:
            rate_limited_task.seconds_per_task = 0

        rate_limited_task.last_call_time = 0
        return rate_limited_task

    return wrap


