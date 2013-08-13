from pyramid.view import view_config
from ptmscout.utils import webutils, decorators
from ptmscout.config import strings
from ptmscout.views.protein import decorators as prot_decorators

@view_config(route_name='protein_summary', renderer='ptmscout:templates/proteins/protein_info.pt')
@prot_decorators.experiment_filter
@decorators.get_protein('id')
def protein_experiment_data_view(context, request, prot):
    return {'pageTitle': strings.protein_summary_page_title % (prot.name),
            'protein': prot}

