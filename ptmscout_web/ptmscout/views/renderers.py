import csv
import StringIO

class TSVRenderer(object):
    def __init__(self, info):
        pass
    
    def __call__(self, value, system):
        """ Call the renderer implementation with the value
        and the system value passed in as arguments and return
        the result (a string or unicode object).  The value is
        the return value of a view.  The system value is a
        dictionary containing available system values
        (e.g. view, context, and request). """
        
        header = value['header']
        data = value['data']
        
        output = StringIO.StringIO()
        writer = csv.writer(output, delimiter = '\t')
        
        writer.writerow(header)
        writer.writerows(data)
        
        rval = output.getvalue()
        output.close()
        
        return rval