from ptmscout.views.protein import structure_view
from pyramid.view import view_config

@view_config(route_name='protein_main', renderer='ptmscout:templates/proteins/protein_viewer.pt')
def protein_view(request):
    return structure_view.protein_structure_viewer(request)

