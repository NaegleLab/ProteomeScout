from pyramid.view import view_config
from ptmscout.database import experiment
from ptmscout.config import strings
from pyramid.httpexceptions import HTTPForbidden

@view_config(route_name='experiment_errors', renderer='ptmscout:/templates/experiments/experiment_errors.pt', permission='private')
def experiment_errors_view(request):
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    user_owner = request.user.experimentOwner(exp)
    
    if not user_owner:
        raise HTTPForbidden()
    
    error_list = [(error.line, error.accession, error.peptide, error.message) for error in exp.errors]
    rejected = len(set( [ error.peptide for error in exp.errors ] ))
    
    return {'user_owner': user_owner,
            'errors':sorted( error_list, key=lambda item: item[0] ),
            'pageTitle': strings.experiment_errors_page_title,
            'experiment':exp,
            'rejected_peptides': rejected,
            'error_count': len(error_list)}