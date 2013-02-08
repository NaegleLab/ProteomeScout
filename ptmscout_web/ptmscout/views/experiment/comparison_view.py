from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications, protein, upload
from ptmscout.utils import forms, webutils, downloadutils, uploadutils
from pyramid.httpexceptions import HTTPFound, HTTPForbidden

@view_config(route_name='experiment_compare', renderer='ptmscout:templates/experiments/experiment_compare.pt')
def experiment_comparison_view(request):
    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, user=request.user)

    return {'pageTitle': strings.experiment_compare_page_title,
            'experiment': exp}
