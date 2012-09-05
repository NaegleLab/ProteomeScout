from pyramid.renderers import get_renderer
from pyramid.events import subscriber
from pyramid.events import BeforeRender

from config import homeUrl, documentationUrl, adminEmail

@subscriber(BeforeRender)
def add_path_definitions(event):
    event['homeUrl'] = homeUrl
    event['documentationUrl'] = documentationUrl
    event['adminEmail'] = adminEmail 
    event['layout'] =  site_layout()

def site_layout():
    renderer = get_renderer("templates/layout.pt")
    layout = renderer.implementation().macros['layout']
    return layout