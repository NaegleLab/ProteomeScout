from types import NoneType

def get(request, var, default):
    return __check_array(request.GET, var, default)
    
def post(request, var, default):
    return __check_array(request.POST, var, default)
    
def __check_array(array, var, default):
    try:
        val = array.getall(var)
        return val[0]
    except KeyError:
        pass
    except IndexError:
        pass
     
    return default