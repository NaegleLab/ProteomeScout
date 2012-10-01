from ptmscout.views.protein.modifications_view import protein_modifications_view
from pyramid.view import view_config


@view_config(route_name='protein_main', renderer='ptmscout:templates/proteins/protein_modifications.pt')
def protein_view(request):
    return protein_modifications_view(request)
