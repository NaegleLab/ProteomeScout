import time
import os

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


