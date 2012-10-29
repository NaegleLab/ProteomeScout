from types import NoneType

def get(request, var, default):
    return __check_array(request.GET, var, default)
    
def post(request, var, default):
    return __check_array(request.POST, var, default)
    
def __check_array(array, var, default):
    try:
        val = array[var]
        return val
    except KeyError:
        pass
     
    return default


def call_catch(errors, method, *args, **kwargs):
    try:
        method(args, kwargs)
    except Exception, e:
        errors.append(e)