from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications
from ptmscout.utils import forms, webutils

def prepopulate_upload_session(exp, user, schema):
    pass

def create_ambiguity_schema(measurements, request):
    pep_list = []
    form_schema = forms.FormSchema()

    i=0
    for ms in measurements:
        fvalues = [ (amb.alt_accession, amb.alt_accession) for amb in ms.ambiguities ]
        field_id = 'ms%d' % (ms.id)
        form_schema.add_select_field(field_id, 'MS %d' % (ms.id), fvalues)
        form_schema.set_required_field(field_id)

        sites = [ (modpep.peptide.getName(), modpep.modification.name) for modpep in ms.peptides ]

        pep_list.append(( i, ms.id, ms.protein.name, ms.protein.species.name, ms.peptide, sites ))
        i+=1

    form_schema.parse_fields(request)

    return form_schema, pep_list

@view_config(route_name='experiment_ambiguity', renderer='ptmscout:templates/experiments/experiment_ambiguity.pt')
def experiment_ambiguity_view(request):
    submitted = webutils.post(request, 'submitted', False) == 'true'

    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, user=request.user)
    peptides = modifications.getMeasuredPeptidesByExperiment(eid, user=request.user)

    form_schema, pep_list = create_ambiguity_schema(peptides, request)
    errors = []

    if submitted:
        errors = forms.FormValidator(form_schema).validate()

        if errors == []:
            prepopulate_upload_session(exp, request.user, form_schema)

    return {'errors': errors,
            'experiment':exp,
            'formrenderer': forms.FormRenderer(form_schema),
            'peptides': pep_list,
            'pageTitle': strings.experiment_ambiguity_page_title}
