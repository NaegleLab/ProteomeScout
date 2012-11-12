from pyramid.view import view_config
from ptmscout.database import experiment
from ptmscout.config import strings


@view_config(route_name='experiment_errors', renderer='ptmscout:/templates/experiments/experiment_errors.pt')
def experiment_errors_view(request):
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    user_owner = exp in request.user.myExperiments()
    error_list = [(error.line, error.accession, error.peptide, error.message) for error in exp.errors]
    rejected = len(set( [ error.peptide for error in exp.errors ] ))
    
    return {'user_owner': user_owner,
            'errors':sorted( error_list, key=lambda item: item[0] ),
            'pageTitle': strings.experiment_errors_page_title,
            'experiment':exp,
            'rejected_peptides': rejected,
            'error_count': len(error_list)}