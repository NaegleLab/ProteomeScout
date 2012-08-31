from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound
import database.experiment as experiment

@view_config(route_name='redirect_to_experiments')
def redirect_to_experiments(request):
    return HTTPFound(request.application_url + "/experiments")

@view_config(route_name='experiments', renderer='templates/experiments.pt')
def experiment_listing(request):
    experiments = experiment.getAllExperiments()
    return {'pageTitle': "Home - Experiments",
            'experiments': experiments}