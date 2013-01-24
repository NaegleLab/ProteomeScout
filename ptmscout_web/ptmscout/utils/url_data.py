import pickle
import base64
import urllib

class URLBuilder(object):
    STRING=0
    INTEGER=1
    OBJECT=2

    def __init__(self):
        self.fields = {}
        self.types = {}
        self.app_url = None

    def create_field(self, field_name, encoding=STRING):
        self.fields[field_name] = None
        self.types[field_name] = encoding

    def decode_request(self, request):
        for field in self.fields:
            if field not in request.GET:
                continue
            value = request.GET[field].strip()

            if self.types[field] == URLBuilder.OBJECT:
                value = pickle.loads( base64.urlsafe_b64decode( value ) )
            elif self.types[field] == URLBuilder.INTEGER:
                value = int(value)
            elif value == "":
                value = None

            self.fields[field] = value

        self.app_url = request.application_url
        self.req_path = request.path

    def __encode_field(self, field):
        value = self.fields[field]
        encoding = self.types[field]

        if value == None:
            return value

        if encoding == URLBuilder.INTEGER:
            value = str(value)
        elif encoding == URLBuilder.OBJECT:
            value = base64.urlsafe_b64encode( pickle.dumps( value ) )

        return value

    def encode_url(self, page):
        if page[0]=='/':
            base_url = self.app_url + page
        else:
            base_url = '/'.join( [self.app_url, page] )

        params = {}
        for field in self.fields.keys():
            v = self.__encode_field(field)
            if v:
                params[field] = v
        if len(params) > 0:
            return '?'.join( [base_url, urllib.urlencode(params)] )

        return base_url

    def get_field(self, field):
        if field not in self.fields:
            raise ValueError("No such field %s" % ( field ) )
        return self.fields[field]


    def field_is_set(self, field):
        return field in self.fields and self.fields[field] != None

    def clear_url(self):
        if self.req_path[0] == '/':
            return self.app_url + self.req_path
        return '/'.join( [self.app_url, self.req_path] )
