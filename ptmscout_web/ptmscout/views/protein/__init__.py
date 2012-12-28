from ptmscout.views.protein.structure_view import protein_structure_viewer
from pyramid.view import view_config


@view_config(route_name='protein_main', renderer='ptmscout:templates/proteins/protein_viewer.pt')
def protein_view(request):
    return protein_structure_viewer(request)
