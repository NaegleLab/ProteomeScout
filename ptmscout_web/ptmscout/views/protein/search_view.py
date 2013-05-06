from ptmscout.config import strings
from ptmscout.database import protein, modifications, taxonomies
from ptmscout.utils import webutils, forms, paginate
from pyramid.view import view_config

QUERY_PAGE_LIMIT=50

def build_schema(request, valid_species):
    submitted = webutils.get(request, 'submitted', False) == 'true'
    form_schema = forms.FormSchema()

    form_schema.add_text_field('acc_search', 'Protein', default='')
    form_schema.add_checkbox_field('include_name', 'Search Protein Names')
    form_schema.add_text_field('pep_search', 'Peptide', default='')
    form_schema.add_text_field('submitted', '', default='false')

    formatted_species = ['all'] + [ sp.capitalize() for sp in valid_species ]
    form_schema.add_autocomplete_field('species', "Species", formatted_species,
            width=100, maxlen=300)
    form_schema.parse_fields(request)

    return submitted, form_schema

def validate_at_least_one_of(items, schema):
    for item in items:
        if schema.field_was_attempted(item):
            return None
    return "At least one of %s are required" % (', '.join([schema.field_names[i] for i in items]))


def build_validator(form_schema):
    validator = forms.FormValidator(form_schema)
    validator.add_alternate_validator('acc_search', lambda name, value, schema: validate_at_least_one_of(['acc_search', 'pep_search'], schema))

    return validator


def perform_query(form_schema, pager, exp_id=None):
    acc_search = form_schema.get_form_value('acc_search')
    include_name = form_schema.get_form_value('include_name') != None
    pep_search = form_schema.get_form_value('pep_search')
    selected_species = form_schema.get_form_value('species')

    if selected_species == 'all' or selected_species == '':
        selected_species=None

    limit, offset = pager.get_pager_limits()
    protein_cnt, proteins = protein.searchProteins(search=acc_search, species=selected_species, sequence=pep_search, page=(limit, offset), exp_id=exp_id, includeNames=include_name)

    pager.set_result_size(protein_cnt)

    return sorted(proteins, key=lambda prot: prot.acc_gene)

def get_protein_metadata(prot, metadata_map, user, exp_id=None):
    measured = modifications.getMeasuredPeptidesByProtein(prot.id, user)

    exp_ids = set()
    residues = set()
    ptms = set()
    sites = set()

    for ms in measured:
        exp_ids.add(ms.experiment_id)

    if exp_id:
        measured = [ ms for ms in measured if ms.experiment_id == exp_id ]

    for ms in measured:
        for mspep in ms.peptides:
            residues.add(mspep.peptide.site_type)
            sites.add(mspep.peptide.site_pos)
            ptms.add(mspep.modification.name)

    metadata_map[prot.id] = ( len(prot.sequence), len(exp_ids), len(sites), ','.join(residues), ', '.join(ptms) )

@view_config(route_name='protein_search', renderer='ptmscout:templates/proteins/protein_search.pt')
def protein_search_view(request):
    species_list = [ species.name for species in taxonomies.getAllSpecies() ]
    submitted, form_schema = build_schema(request, species_list)

    pager = paginate.Paginator(form_schema, QUERY_PAGE_LIMIT)
    pager.parse_parameters(request)

    proteins = []
    protein_metadata = {}
    errors = []

    if submitted:
        errors = build_validator(form_schema).validate()
        if len(errors) == 0:
            proteins = perform_query(form_schema, pager)

    for p in proteins:
        get_protein_metadata(p, protein_metadata, request.user)

    form_renderer = forms.FormRenderer(form_schema)
    return {'pageTitle': strings.protein_search_page_title,
            'formrenderer':form_renderer,
            'paginator': pager,
            'proteins':proteins,
            'protein_metadata':protein_metadata,
            'errors': errors,
            'submitted': submitted,
            'search_url': "%s/proteins" % request.application_url}
