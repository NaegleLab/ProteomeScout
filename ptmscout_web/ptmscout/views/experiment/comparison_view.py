from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications, protein, upload
from ptmscout.utils import forms, webutils, downloadutils, uploadutils
from pyramid.httpexceptions import HTTPFound, HTTPForbidden

def format_peptide_list(peptide_list):
    formatted_peptides = []
    for ms, modpep in peptide_list:
        formatted = {'id':ms.id,
                     'gene':ms.protein.getGeneName(),
                     'protein':ms.protein.name, 
                     'tryps':ms.peptide,
                     'align':modpep.peptide.pep_aligned,
                     'site':modpep.peptide.getName(),
                     'mod':modpep.modification.name }

        formatted_peptides.append(formatted)
    return formatted_peptides

def compare_to_all(exp, user, experiment_list = set()):
    exps = experiment.getAllExperiments(user)
    if experiment_list == set():
        experiment_list = set([e.id for e in exps])

    experiment_info = {}
    for exp in exps:
        experiment_info[exp.id] = { 'name': exp.name, 'export': exp.export }

    by_experiment = {}
    for eid in experiment_list:
        by_experiment[eid] = []

    ambiguous_peptides = []
    novel_sites = []

    for ms in exp.measurements:
        if ms.isAmbiguous():
            for modpep in ms.peptides:
                ambiguous_peptides.append( (ms, modpep) )
            continue

        for modpep in ms.peptides:
            other_exps = modifications.getExperimentsReportingModifiedPeptide(modpep, exps)
            if len(other_exps) == 1 and other_exps[0].id == exp.id:
                novel_sites.append( (ms, modpep) )

            for exp in other_exps:
                by_experiment[exp.id].append( (ms, modpep) )

    for exp_id in experiment_list:
        by_experiment[exp_id] = format_peptide_list( by_experiment[exp_id] )

    return {'ambiguous': format_peptide_list(ambiguous_peptides), 'novel': format_peptide_list(novel_sites), 'by_experiment': by_experiment, 'experiment_info': experiment_info}

@view_config(route_name='experiment_compare', renderer='ptmscout:templates/experiments/experiment_compare.pt')
def experiment_comparison_view(request):
    submitted_val = webutils.post(request, 'submitted', False)

    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, user=request.user)

    results = None
    if submitted_val == 'all':
        results = compare_to_all(exp, request.user)
    elif submitted_val == 'subset':
        experiment_list = set([int(eid) for eid in request.POST.getall('experiment')])
        results = compare_to_all(exp, request.user, experiment_list)


    return {'pageTitle': strings.experiment_compare_page_title,
            'experiment': exp,
            'results': results}
