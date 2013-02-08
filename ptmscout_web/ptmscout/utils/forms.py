import cgi
import webutils
from ptmscout.config import strings

class FormLiteral(object):
    def __init__(self, html):
        self.html = html
    
    def __html__(self):
        return self.html

field_not_empty_test = lambda field_value: field_value != None and field_value != ''


def field_equals_test(value):
    return lambda field_value: field_value == value

class FormSchema(object):
    CHECKBOX=1
    TEXT=2
    TEXTAREA=3
    SELECT=4
    RADIO=5
    PASSWORD=6
    FILE=7
    
    def __init__(self):
        self.form_values = {}
        
        self.field_names = {}
        self.field_opts = {}
        self.field_defaults = {}
        self.field_types = {}
        
        self.enum_fields = set()
        self.enum_values = {}
        
        self.numeric_fields = set()
        self.required_fields = set()
        self.conditional_fields = {}

    def set_field(self, ref, value):
        self.form_values[ref] = value

    def get_form_value(self, ref):
        value = self.form_values[ref]
        if value == None:
            value = self.field_defaults[ref]
        return value
    
    def field_was_attempted(self, field_ref):
        field_value = self.form_values[field_ref]
        return not (field_value == None or field_value == '')

    def parse_fields(self, request):
        for ref in self.field_names:
            v = webutils.post(request, ref, None) 
            if v == None: v = webutils.get(request, ref, None)
            
            if v != None and self.field_types[ref] != FormSchema.FILE:
                v = v.strip()
            self.form_values[ref] = v

    def parse_fields_dict(self, field_dict):
        for ref in self.field_names:
            v = None
            if ref in field_dict:
                v = field_dict[ref]
            if v != None and self.field_types[ref] != FormSchema.FILE:
                v = str(v)
                v = v.strip()
            self.form_values[ref] = v
   
    def is_defaulted(self):
        all_defaults = [ self.form_values[ref] == self.field_defaults[ref] for ref in self.form_values ]
        return reduce(bool.__and__, all_defaults, True)

    def set_field_required_condition(self, ref, parent, condition):
        self.conditional_fields[ref] = (parent, condition)
    
    def set_required_field(self, ref):
        self.required_fields.add(ref)
    
    
    def add_checkbox_field(self, ref, name, default=None):
        self.field_names[ref] = name
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.CHECKBOX
    
    def add_numeric_field(self, ref, name, maxlen=None, default=None):
        self.field_names[ref] = name
        self.numeric_fields.add(ref)
        
        width=None if maxlen==None else maxlen+1
        self.field_opts[ref] = (maxlen, width,)
        
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.TEXT
    
    def add_text_field(self, ref, name, maxlen=None, width=None, default=None):
        self.field_names[ref] = name
        self.field_opts[ref] = (maxlen, width,)
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.TEXT
        
    def add_password_field(self, ref, name, maxlen=None, width=None, default=None):
        self.field_names[ref] = name
        self.field_opts[ref] = (maxlen, width,)
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.PASSWORD

    def add_textarea_field(self, ref, name, width, height, default=None):
        self.field_names[ref] = name
        self.field_opts[ref] = (width, height)
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.TEXTAREA
    
    def add_radio_field(self, ref, name, fvalues=[], default=None):
        self.field_names[ref] = name
        self.enum_values[ref] = fvalues
        self.enum_fields.add(ref)
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.RADIO 
    
    def add_select_field(self, ref, name, fvalues=[], default=None):
        self.field_names[ref] = name
        self.enum_values[ref] = fvalues
        self.enum_fields.add(ref)
        self.field_defaults[ref] = default
        self.field_types[ref] = FormSchema.SELECT

    def add_file_upload_field(self, ref, name):
        self.field_names[ref] = name
        self.field_defaults[ref] = ''
        self.field_types[ref] = FormSchema.FILE


class FormRenderer(object):
    def __init__(self, schema):
        self.schema = schema
        self.radio_indexes = {}
    
    def render(self, ref, id_=None, class_=None):
        if ref not in self.schema.field_names:
            raise Exception("No such form field: %s" % ref)
        
        field_type = self.schema.field_types[ref]
        id_str = '' if id_ == None else 'id="%s"' % (id_)
        cls_str = '' if class_ == None else 'class="%s"' % (class_)
        
        if field_type == FormSchema.TEXT:
            return self.__render_text(ref, id_str, cls_str)
        if field_type == FormSchema.TEXTAREA:
            return self.__render_textarea(ref, id_str, cls_str)
        if field_type == FormSchema.SELECT:
            return self.__render_select(ref, id_str, cls_str)
        if field_type == FormSchema.CHECKBOX:
            return self.__render_checkbox(ref, id_str, cls_str)
        if field_type == FormSchema.RADIO:
            return self.__render_radio(ref, id_str, cls_str)
        if field_type == FormSchema.PASSWORD:
            return self.__render_password(ref, id_str, cls_str)
        if field_type == FormSchema.FILE:
            return self.__render_file(ref, id_str, cls_str)
    
    def __render_select(self, ref, id_str, cls_str):
        items = []
        items.append('<select %s %s name="%s">' % (id_str, cls_str, ref))
        
        cur_value = self.schema.get_form_value(ref)
        if cur_value == None:
            cur_value = ''

        valid_values = self.schema.enum_values[ref][:]
        if self.schema.field_defaults[ref] == None:
            valid_values = [('','')] + valid_values

        for (value, proper_name) in valid_values:
            selected = 'selected' if cur_value == value else ''
            items.append('<option value="%s" %s>%s</option>' % (value, selected, proper_name))
        
        items.append('</select>')
        
        return FormLiteral("\n".join(items))
    
    def __render_checkbox(self, ref, id_str, cls_str):
        checked = 'checked' if self.schema.form_values[ref] != None else '' 
        proper_name = self.schema.field_names[ref]
        
        return FormLiteral('<input type="checkbox" %s %s name="%s" %s /> %s' % (id_str, cls_str, ref, checked, proper_name))
    
    def __render_radio(self, ref, id_str, cls_str):
        if ref not in self.radio_indexes:
            self.radio_indexes[ref] = 0
        
        index = self.radio_indexes[ref]
        cur_value = self.schema.get_form_value(ref)
        
        if index < len(self.schema.enum_values[ref]):
            (value, proper_name) = self.schema.enum_values[ref][index]
            selected = 'checked' if cur_value == value else ''
            self.radio_indexes[ref] = index+1
            return FormLiteral('<input type="radio" %s %s name="%s" value="%s" %s /> %s' % (id_str, cls_str, ref, value, selected, proper_name))
        
        return FormLiteral('')
    
    def __render_text(self, ref, id_str, cls_str):
        value = self.schema.get_form_value(ref)
        value_str = '' if value == None else 'value="%s"' % (cgi.escape(value),)
        
        maxlen, size = self.schema.field_opts[ref]
        
        size_str = '' if size == None else 'size="%d"' % (size)
        maxlen_str = '' if maxlen == None else 'maxlength="%d"' % (maxlen)
        
        html = '<input type="text" %s %s %s %s name="%s" %s />' % (id_str, cls_str, size_str, maxlen_str, ref, value_str)
        return FormLiteral(html)
    
    def __render_textarea(self, ref, id_str, cls_str):
        value = self.schema.get_form_value(ref)
        if value == None:
            value = ''
            
        width, height = self.schema.field_opts[ref]
        
        html = '<textarea cols="%d" rows="%d" name="%s">%s</textarea>' % (width, height, ref, cgi.escape(value))
        return FormLiteral(html)
    
    def __render_password(self, ref, id_str, cls_str):
        maxlen, size = self.schema.field_opts[ref]
        
        size_str = '' if size == None else 'size="%d"' % (size)
        maxlen_str = '' if maxlen == None else 'maxlength="%d"' % (maxlen)
        
        html = '<input type="password" %s %s %s %s name="%s" />' % (id_str, cls_str, size_str, maxlen_str, ref)
        return FormLiteral(html)
    
    def __render_file(self, ref, id_str, cls_str):
        html='<input %s %s type="file" name="%s" />' % (id_str, cls_str, ref)
        return FormLiteral(html)


class FormValidator(object):
    def __init__(self, schema):
        self.schema = schema
        self.extra_validators = {}

    # validators are functions or lambda expressions accepting 3 args: field_name, field_value, schema
    def add_alternate_validator(self, ref, validator):
        self.extra_validators[ref] = validator

    def validate(self):
        errors = []
        
        for ref in self.schema.field_names:
            value = self.schema.get_form_value(ref)
            error = self.validate_field(ref, value)
            
            if error != None:
                errors.append(error)
        
        return errors
    
    def validate_field(self, field_ref, field_value):
        proper_name = self.schema.field_names[field_ref]
        
        if field_ref in self.schema.required_fields and not self.schema.field_was_attempted(field_ref):
            return strings.failure_reason_required_fields_cannot_be_empty % (proper_name)

        if field_ref in self.schema.conditional_fields:
            parent_ref, condition = self.schema.conditional_fields[field_ref]
            parent_value = self.schema.get_form_value(parent_ref)
            
            if condition(parent_value) and not self.schema.field_was_attempted(field_ref):
                return strings.failure_reason_required_fields_cannot_be_empty % (proper_name)
        
        if field_ref in self.schema.numeric_fields and self.schema.field_was_attempted(field_ref):
            try:
                int(field_value)
            except:
                return strings.failure_reason_field_must_be_numeric % (proper_name)
        
        if field_ref in self.schema.enum_fields and self.schema.field_was_attempted(field_ref):
            internal_values = set([ k for k, _ in self.schema.enum_values[field_ref] ])
            if field_value not in internal_values:
                return strings.failure_reason_field_value_not_valid % (proper_name)
            
        if field_ref in self.extra_validators:
            validator = self.extra_validators[field_ref]
            error = validator(proper_name, field_value, self.schema)
            if error != None:
                return error
            
        return None



