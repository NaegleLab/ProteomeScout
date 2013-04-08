import json, base64

def generate_metadata_query(field, value):
    query = [ ['nop', ['metadata', field, 'eq', value] ] ]
    return base64.b64encode( json.dumps( query ) )

def generate_scansite_query(category, value):
    query = [ ['nop', ['scansite', category, 'eq', value, '5'] ] ]
    return base64.b64encode( json.dumps( query ) )
