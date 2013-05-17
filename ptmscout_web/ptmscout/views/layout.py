from pyramid.renderers import get_renderer
from pyramid.events import subscriber
from pyramid.events import BeforeRender

from ptmscout.config.settings import documentationUrl, adminEmail, naegle_lab_url

@subscriber(BeforeRender)
def add_path_definitions(event):
    event['documentationUrl'] = documentationUrl
    event['naegleUrl'] = naegle_lab_url
    event['adminEmail'] = adminEmail 
    event['parent_link'] = None
    event['layout'] = site_layout()
    event['experiment'] = None

    event['experiment_header'] = experiment_template()
    event['experiment_info'] = experiment_info()
    event['experiment_summary_menu'] = experiment_summary_template()    
    event['protein_header'] = protein_template()
    event['protein_list'] = protein_list_template()
    event['protein_info'] = protein_info_template()
    event['protein_nav'] = protein_nav_template()
    event['protein_search_form'] = protein_search_form()
    
    event['dataset_explorer'] = dataset_explorer_template()

    event['footer'] = footer_template()
    
    event['redirect'] = None

def site_layout():
    renderer = get_renderer("ptmscout:templates/macro/layout.pt")
    layout = renderer.implementation().macros['layout']
    return layout

def experiment_template():
    renderer = get_renderer("ptmscout:templates/macro/experiment_header.pt")
    legend = renderer.implementation().macros['experiment_header']
    return legend

def experiment_info():
    renderer = get_renderer("ptmscout:templates/macro/experiment_info.pt")
    legend = renderer.implementation().macros['experiment_info']
    return legend

def experiment_summary_template():
    renderer = get_renderer("ptmscout:templates/macro/experiment_summary_menu.pt")
    legend = renderer.implementation().macros['experiment_summary_menu']
    return legend

def protein_template():
    renderer = get_renderer("ptmscout:templates/macro/protein_header.pt")
    legend = renderer.implementation().macros['protein_header']
    return legend

def protein_list_template():
    renderer = get_renderer("ptmscout:templates/macro/protein_list.pt")
    legend = renderer.implementation().macros['protein_list']
    return legend

def protein_info_template():
    renderer = get_renderer("ptmscout:templates/macro/protein_info.pt")
    legend = renderer.implementation().macros['protein_info']
    return legend

def protein_search_form():
    renderer = get_renderer("ptmscout:templates/macro/protein_search_form.pt")
    legend = renderer.implementation().macros['protein_search_form']
    return legend

def protein_nav_template():
    renderer = get_renderer("ptmscout:templates/macro/protein_nav.pt")
    legend = renderer.implementation().macros['protein_nav']
    return legend

def dataset_explorer_template():
    renderer = get_renderer("ptmscout:templates/dataset/dataset_explorer.pt")
    legend = renderer.implementation().macros['dataset_explorer']
    return legend

def footer_template():
    renderer = get_renderer("ptmscout:templates/macro/footer.pt")
    legend = renderer.implementation().macros['footer']
    return legend
