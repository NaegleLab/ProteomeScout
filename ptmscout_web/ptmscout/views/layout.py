from pyramid.renderers import get_renderer
from pyramid.events import subscriber
from pyramid.events import BeforeRender

from ptmscout.config.settings import homeUrl, documentationUrl, adminEmail

@subscriber(BeforeRender)
def add_path_definitions(event):
    event['homeUrl'] = homeUrl
    event['documentationUrl'] = documentationUrl
    event['adminEmail'] = adminEmail 
    event['parent_link'] = None
    event['layout'] = site_layout()
    event['experiment'] = None
    event['experiment_header'] = experiment_template()
    event['protein_header'] = protein_template()
    event['protein_list'] = protein_list_template()
    event['redirect'] = None

def site_layout():
    renderer = get_renderer("ptmscout:templates/macro/layout.pt")
    layout = renderer.implementation().macros['layout']
    return layout

def experiment_template():
    renderer = get_renderer("ptmscout:templates/macro/experiment_header.pt")
    legend = renderer.implementation().macros['experiment_header']
    return legend

def protein_template():
    renderer = get_renderer("ptmscout:templates/macro/protein_header.pt")
    legend = renderer.implementation().macros['protein_header']
    return legend

def protein_list_template():
    renderer = get_renderer("ptmscout:templates/macro/protein_list.pt")
    legend = renderer.implementation().macros['protein_list']
    return legend