from pyramid.view import view_config
from ptmscout.config import strings
from ptmscout.database import experiment, modifications, protein, upload
from ptmscout.utils import forms, webutils, downloadutils, uploadutils
from pyramid.httpexceptions import HTTPFound, HTTPForbidden
import re

def create_session(exp, user, exp_file, name_prefix, columns, units):
    session = upload.Session()
    
    session.user_id = user.id
    session.data_file = exp_file
    session.load_type = 'extension'
    session.stage = 'metadata'
    session.change_description = ''
    
    new_exp = experiment.Experiment()
    new_exp.copyData(exp)
    new_exp.name = name_prefix + new_exp.name
    new_exp.saveExperiment()

    session.experiment_id = new_exp.id
    session.parent_experiment = exp.id
    session.change_description = "New peptide to protein accession assignments"

    numcols = len(columns)
    session.columns = []
    for c in xrange(0, numcols):
        col = upload.SessionColumn()
        col.type = columns[c]['type']
        col.label = columns[c]['label']
        col.column_number = c
        session.columns.append(col)
    session.units = units

    session.save()
    
    return session.id

def prepopulate_upload_session(exp, user, schema, was_defaults):
    ms_map = {}
    for field in schema.field_names:
        m = re.match("ms([0-9]+)", field)
        ms_id = int(m.group(1))
        ms_map[ms_id] = schema.get_form_value(field)

    headers, exp_filename = downloadutils.experiment_to_tsv(exp, ms_map)
    assignments = uploadutils.assign_columns_by_name(headers)

    if was_defaults:
        name_prefix = '[Default Assignments] '
    else:
        name_prefix = '[Changed Assignments] '

    session_id = create_session(exp, user, exp_filename, name_prefix, assignments['columns'], assignments['units'])
    return session_id

def assign_defaults(measurements, schema, user):
    changed_msids = set()
    for ms in measurements:
        found_current = False
        protein_choices = []
        for amb in ms.ambiguities:
            alt_proteins = protein.getProteinsByAccession(amb.alt_accession)
            max_prot = None
            max_mods = 0

            for p in alt_proteins:
                mods = len(modifications.getMeasuredPeptidesByProtein(p.id, user))
                if mods > max_mods:
                    max_prot = p
                    max_mods = mods

            if max_prot:
                go_terms = len(max_prot.GO_terms)
                if max_prot.id == ms.protein.id:
                    found_current=True
                protein_choices.append((max_prot, amb.alt_accession, max_mods, go_terms))

        if not found_current:
            mods = len(modifications.getMeasuredPeptidesByProtein(ms.protein_id, user))
            protein_choices.append((ms.protein, ms.query_accession, mods, len(ms.protein.GO_terms)))

        protein_choices = sorted(protein_choices, key=lambda item: (item[2], item[3]))

        chosen = protein_choices[-1][1]

        if chosen != ms.query_accession:
            changed_msids.add(ms.id)

        schema.set_field('ms%d'%(ms.id), chosen)

    return changed_msids

def create_ambiguity_schema(measurements, request):
    pep_list = []
    form_schema = forms.FormSchema()

    i=0
    for ms in measurements:
        alt_accessions = set([amb.alt_accession for amb in ms.ambiguities])
        fvalues = [ (amb.alt_accession, "%s : %s" % (amb.locus, amb.alt_accession)) for amb in ms.ambiguities ]
        if ms.query_accession not in alt_accessions:
            fvalues = [ (ms.query_accession, "%s : %s" % (ms.protein.acc_gene, ms.query_accession)) ] + fvalues

        field_id = 'ms%d' % (ms.id)
        form_schema.add_select_field(field_id, 'MS %d' % (ms.id), fvalues, default=ms.query_accession)
        form_schema.set_required_field(field_id)

        sites = [ (modpep.peptide.getName(), modpep.modification.name) for modpep in ms.peptides ]

        pep_list.append(( i, ms.id, ms.query_accession, ms.protein.acc_gene, ms.protein.name, ms.protein.species.name, ms.peptide, sites ))
        i+=1

    form_schema.parse_fields(request)

    return form_schema, pep_list

@view_config(route_name='experiment_ambiguity', renderer='ptmscout:templates/experiments/experiment_ambiguity.pt')
def experiment_ambiguity_view(request):
    submitted = webutils.post(request, 'submitted', False) == 'true'
    defaults = webutils.get(request, 'defaults', False) == 'true'

    eid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(eid, user=request.user)

    if exp.ambiguity == 0:
        raise HTTPForbidden()

    peptides = modifications.getMeasuredPeptidesByExperiment(eid, user=request.user)

    form_schema, pep_list = create_ambiguity_schema(peptides, request)

    changed_msids = set()
    if defaults and not submitted:
        changed_msids = assign_defaults(peptides, form_schema, request.user)

    errors = []

    if submitted:
        errors = forms.FormValidator(form_schema).validate()

        if form_schema.is_defaulted():
            errors.append(strings.experiment_ambiguity_error_no_change)

        if errors == []:
            session_id = prepopulate_upload_session(exp, request.user, form_schema, defaults)
            raise HTTPFound('%s/upload/%d' % (request.application_url, session_id))

    return {'errors': errors,
            'experiment':exp,
            'formrenderer': forms.FormRenderer(form_schema),
            'peptides': pep_list,
            'assigned_defaults': defaults,
            'changed_default': changed_msids,
            'pageTitle': "%s: %s" % (strings.experiment_ambiguity_page_title, exp.name)}
