from pyramid.renderers import get_renderer
from pyramid.events import subscriber
from pyramid.events import BeforeRender

from config import homeUrl, documentationUrl, adminEmail

@subscriber(BeforeRender)
def add_path_definitions(event):
    event['homeUrl'] = homeUrl
    event['documentationUrl'] = documentationUrl
    event['adminEmail'] = adminEmail 
    event['parent_link'] = None
    event['layout'] = site_layout()
    event['experiment_header'] = experiment_template()
    event['redirect'] = None

def site_layout():
    renderer = get_renderer("templates/layout.pt")
    layout = renderer.implementation().macros['layout']
    return layout

def experiment_template():
    renderer = get_renderer("templates/experiment_header.pt")
    legend = renderer.implementation().macros['experiment_header']
    return legend
