from pyramid.view import view_config
import json
import base64
from ptmscout.database import experiment
from ptmscout.utils import protein_utils, motif

def get_pfam_metadata(measurements, metadata_fields):
    metadata_fields.update({'Pfam-Domain':set(),'Pfam-Site':set(),'Region':set()})
    
    for ms in measurements:
        for domain in ms.protein.domains:
            metadata_fields['Pfam-Domain'].add(domain.label)
            
        for modpep in ms.peptides:
            if modpep.peptide.protein_domain:
                metadata_fields['Pfam-Site'].add(modpep.peptide.protein_domain.label)
            
        for region in ms.protein.regions:
            for modpep in ms.peptides:
                if region.hasSite(modpep.peptide.site_pos):
                    metadata_fields['Region'].add(region.label)

def get_protein_metadata(measurements, metadata_fields, accessions):
    metadata_fields = {
                       'Gene':set(),
                       'Protein Accession':set(),
                       'Protein Name':set(),
                       'Species':set(),
                       }
    
    for ms in measurements:
        metadata_fields['Gene'].add(ms.protein.getGeneName())
        metadata_fields['Protein Name'].add(ms.protein.name)
        metadata_fields['Species'].add(ms.protein.species.name)
        
        for acc in ms.protein.accessions:
            metadata_fields['Protein Accession'].add(acc.value)
        
    accessions += sorted( list(metadata_fields['Gene'] | metadata_fields['Protein Name'] | metadata_fields['Protein Accession']) )
    del metadata_fields['Gene']
    del metadata_fields['Protein Name']
    del metadata_fields['Protein Accession']
        
def get_scansite_metadata(measurements, metadata_fields):
    metadata_fields.update({'Scansite-Kinase':set(),
                            'Scansite-Bind':set()})
    
    for ms in measurements:
        for modpep in ms.peptides:
            for scansite in modpep.peptide.predictions:
                if scansite.source == 'scansite_kinase':
                    metadata_fields['Scansite-Kinase'].add(scansite.value)
                if scansite.source == 'scansite_bind':
                    metadata_fields['Scansite-Bind'].add(scansite.value)

    for key in metadata_fields:
        metadata_fields[key] = list(metadata_fields[key])


def get_GO_metadata(measurements, metadata_fields):
    metadata_fields.update({'GO-Cellular Component':set(),
                            'GO-Molecular Function':set(),
                            'GO-Biological Process':set()})

    for ms in measurements:
        for go_term in ms.protein.GO_terms:
            if go_term.GO_term.aspect == 'C':
                metadata_fields['GO-Cellular Component'].add(go_term.GO_term.fullName())
            if go_term.GO_term.aspect == 'P':
                metadata_fields['GO-Biological Process'].add(go_term.GO_term.fullName())
            if go_term.GO_term.aspect == 'F':
                metadata_fields['GO-Molecular Function'].add(go_term.GO_term.fullName())


def get_modification_metadata(measurements, metadata_fields):
    metadata_fields.update({'Modified Residue':set(),
                            'Modification Type':set()})

    for ms in measurements:
        for modpep in ms.peptides:
            metadata_fields['Modified Residue'].add(modpep.peptide.site_type)
            
            metadata_fields['Modification Type'].add(modpep.modification.name)
            for ptm in modpep.modification.getAllParents():
                metadata_fields['Modification Type'].add(ptm.name)


def format_explorer_view(experiment_id, measurements):
    quantitative_fields = set()
    
    accessions = []
    metadata_fields = {}
    scansite_fields = {}
    cluster_labels = {}
    subset_labels = []
    
    get_protein_metadata(measurements, metadata_fields, accessions)
    get_pfam_metadata(measurements, metadata_fields)
    
    get_scansite_metadata(measurements, scansite_fields)
    
    get_GO_metadata(measurements, metadata_fields)
    get_modification_metadata(measurements, metadata_fields)
    
    for ms in measurements:
        for d in ms.data:
            quantitative_fields.add(d.formatted_label)

    
    for field in metadata_fields.keys():
        metadata_fields[field] = sorted( list(metadata_fields[field]) )
    
    
    field_data = {
                  'experiment_id': experiment_id,
                  'quantitative_fields': sorted( list(quantitative_fields) ),
                  'accessions': accessions,
                  'scansite_fields': scansite_fields,
                  'scansite_keys': sorted( scansite_fields.keys() ),
                  'metadata_fields': metadata_fields,
                  'metadata_keys': sorted( metadata_fields.keys() ),
                  'clustering_sets': sorted( list(cluster_labels.keys()) ),
                  'clustering_labels': cluster_labels,
                  'subset_labels':subset_labels}
    
    return {'field_data': base64.b64encode( json.dumps(field_data) )}

class QueryError(Exception):
    pass

def calculate_feature_enrichment(foreground, background):
    pass


def create_protein_filter(filter_value):
    def protein_filter(ms):
        return ms.protein.hasAccession(filter_value) or \
                ms.protein.name == filter_value or \
                ms.protein.acc_gene == filter_value or \
                ms.protein.locus == filter_value
    return protein_filter


def create_metadata_filter(field, op, value):
    def op_apply_wrapper(filter_func):
        def op_apply(ms):
            result = reduce(bool.__or__, filter_func(ms), False)
            return (op == 'eq' and result) or \
                    (op == 'neq' and not result)
                    
        return op_apply
    
    @op_apply_wrapper
    def modtype_filter(ms):
        return [ modpep.modification.name.lower() == value.lower() for modpep in ms.peptides ]
    
    @op_apply_wrapper    
    def modsite_filter(ms):
        return [ modpep.peptide.site_type == value for modpep in ms.peptides ]
    
    @op_apply_wrapper
    def GO_filter(ms):
        return [ value == go_term.GO_term.fullName() for go_term in ms.protein.GO_terms ]

    @op_apply_wrapper
    def pfam_domain_filter(ms):
        return [ value.lower() == domain.label.lower() for domain in ms.protein.domains ]
    
    @op_apply_wrapper
    def pfam_site_filter(ms):
        return [ None != modpep.peptide.protein_domain.lower() and modpep.peptide.protein_domain.label.lower() == value.lower() for modpep in ms.peptides ]
    
    @op_apply_wrapper
    def region_filter(ms):
        return [ region.label.lower() == value.lower() and region.hasSite(modpep.peptide.site_pos) for region in ms.protein.regions for modpep in ms.peptides ]
                
    if field == 'Modification Type':
        return modtype_filter
    if field == 'Modified Residue':
        return modsite_filter
    if field.find('GO-') == 0:
        return GO_filter
    if field == 'Pfam-Domain':
        return pfam_domain_filter
    if field == 'Pfam-Site':
        return pfam_site_filter
    if field == 'Region':
        return region_filter
    
    raise QueryError("Unrecognized metadata filter field '%s'" % (field))

def create_scansite_filter(field, op, value, stringency):
    try:
        stringency = float(stringency)
    except:
        raise QueryError("Stringency filter must be numeric value")
    
    def op_apply_wrapper(filter_func):
        def op_apply(ms):
            result = reduce(bool.__or__, filter_func(ms), False)
            return (op == 'eq' and result) or \
                    (op == 'neq' and not result)
        return op_apply
    
    @op_apply_wrapper
    def kinase_filter(ms):
        return [prediction.source == 'scansite_kinase' and prediction.value == value and prediction.score < stringency 
                 for modpep in ms.peptides for prediction in modpep.peptide.predictions]
    
    @op_apply_wrapper
    def bind_filter(ms):
        return [prediction.source == 'scansite_bind' and prediction.value == value and prediction.score < stringency 
                 for modpep in ms.peptides for prediction in modpep.peptide.predictions]

    
    if field=='Scansite-Kinase':
        return kinase_filter
    if field=='Scansite-Bind':
        return bind_filter
    
    raise QueryError("Unrecognized scansite filter field '%s'" % (field))

def create_sequence_filter(sequence):
    def seq_filter(ms):
        return reduce(bool.__or__, [modpep.peptide.pep_aligned.upper() == sequence.upper() for modpep in ms.peptides ])

    return seq_filter

def find_in(array, checkset):
        return [i for i, item in enumerate(array) if item in checkset]

def readable_expression(expression):
    readable_op_map = {'eq':"=", 'neq':"\u2260", 'gt':">", 'geq':"\u2265", 'lt':"<", 'leq':"\u2264", 'add':"+", 'sub':"-", 'mul':"x", 'div':"/"}
    return [ readable_op_map[term] if term in readable_op_map else term 
                for term in expression ]

class MissingMeasurementError(Exception):
    pass

def parse_value(value):
    try:
        num = float(value)
        def get_data(ms):
                return num
        return get_data
    except:
        def get_data(ms):
            for d in ms.data:
                if d.formatted_label == value:
                    if d.value == None:
                        raise MissingMeasurementError()
                    return d.value
            
            raise QueryError("No such data label: '%s'" % (value))
        return get_data


def parse_arthmetic_expression(expression):
    if len(expression) == 1:
        return parse_value(expression[0])
    elif len(expression) == 3:
        v1 = parse_value(expression[0])
        op = expression[1]
        v2 = parse_value(expression[2])
        
        def adder(ms):
            return v1(ms) + v2(ms)
        def suber(ms):
            return v1(ms) - v2(ms)
        def muler(ms):
            return v1(ms) * v2(ms)
        def diver(ms):
            return v1(ms) / v2(ms)
        
        if op == 'add':
            return adder
        if op == 'sub':
            return suber
        if op == 'mul':
            return muler
        if op == 'div':
            return diver
        
        raise QueryError("Invalid operator type: '%s'" % (op))
        

def create_quantitative_filter(expression):
    comparison_operators = set([ 'eq', 'neq', 'gt', 'lt', 'geq', 'leq' ])
    
    comp_index = find_in(expression, comparison_operators)
    if len(comp_index) == 0:
        raise QueryError("No comparison found in quantitative expression: %s" % (readable_expression(expression)))
    comp_index = comp_index[0]
    
    LHS = expression[:comp_index]
    RHS = expression[comp_index+1:]
    
    LHS_func = parse_arthmetic_expression(LHS)
    RHS_func = parse_arthmetic_expression(RHS)
    
    def catch_wrapper(func):
        def catcher(ms):
            try:
                return func(ms)
            except MissingMeasurementError:
                return False
        return catcher
    
    @catch_wrapper
    def eq_func(ms):
        return LHS_func(ms) == RHS_func(ms)
    
    @catch_wrapper
    def neq_func(ms):
        return LHS_func(ms) != RHS_func(ms)
    
    @catch_wrapper
    def gt_func(ms):
        return LHS_func(ms) > RHS_func(ms)
    
    @catch_wrapper
    def lt_func(ms):
        return LHS_func(ms) < RHS_func(ms)
    
    @catch_wrapper
    def geq_func(ms):
        return LHS_func(ms) >= RHS_func(ms)
    
    @catch_wrapper
    def leq_func(ms):
        return LHS_func(ms) <= RHS_func(ms)
    
    if expression[comp_index] == 'eq':
        return eq_func
    if expression[comp_index] == 'neq':
        return neq_func
    if expression[comp_index] == 'gt':
        return gt_func
    if expression[comp_index] == 'lt':
        return lt_func
    if expression[comp_index] == 'geq':
        return geq_func
    if expression[comp_index] == 'leq':
        return leq_func
    
    raise QueryError("Invalid comparison operator '%s'" % ( expression[comp_index] ))

def parse_condition(condition):
    if condition[0] == 'protein':
        return create_protein_filter(condition[1])
    if condition[0] == 'metadata':
        return create_metadata_filter(*tuple(condition[1:]))
    if condition[0] == 'scansite':
        return create_scansite_filter(*tuple(condition[1:]))
    if condition[0] == 'sequence':
        return create_sequence_filter(condition[1])
    if condition[0] == 'quantitative':
        return create_quantitative_filter(condition[1:])


def parse_dnf_clause(conditions):
    def conjunction(operands):
        def apply_conjunction(ms):
            for op in operands:
                if not op(ms):
                    return False
            return True
        
        return apply_conjunction
    
    return conjunction( [ parse_condition(cond) for cond in conditions] )

def parse_query_expression(expression):
    def passthrough(ms):
        return True
    
    def disjunction(operands):
        def apply_disjunction(ms):
            for op in operands:
                if op(ms):
                    return True
            return False
        
        return apply_disjunction
    
    if expression == 'experiment':
        return passthrough
    
    dnf_clauses = []
    
    for op, condition in expression:
        if op != 'and':
            dnf_clauses.append([])
        
        dnf_clauses[-1].append(condition)
        
    return disjunction( [ parse_dnf_clause(op) for op in dnf_clauses ] )

def format_peptide_data(measurements):
    formatted_peptide_data = []
    for ms in measurements:
        gene_name = ms.protein.getGeneName()
        
        for modpep in ms.peptides:
            pep_aligned = modpep.peptide.pep_aligned
            site_name = modpep.peptide.getName()
            mod_name = modpep.modification.name
            
            formatted_peptide_data.append( (ms.id, ms.query_accession, gene_name, pep_aligned, site_name, mod_name) )
    
    return formatted_peptide_data

def format_measurement_data(measurements):
    experiment_data = []
    
    for ms in measurements:
        series_label = '%s:%s' % (ms.protein.getGeneName(), ','.join([ modpep.peptide.getName() for modpep in ms.peptides ]))
        run_labels = {}
        run_units = {}
        
        for d in ms.data:
            datapts = run_labels.get(d.run, [])
            run_units[d.run] = d.units
            
            datapts.append( (d.priority, d.label, d.value, d.type) )
            run_labels[d.run] = datapts
        
        for run in run_labels.keys():
            experiment_data.append({'name':run, 'units': run_units[run], 'label':series_label, 'series': [(l,v,t) for _,l,v,t in  sorted( run_labels[run] )]})
        
    return experiment_data


@view_config(route_name='compute_subset', renderer='json', permission='private')
def compute_subset_enrichment(request):
    query_expression = json.loads(request.body, encoding=request.charset)
    exp_id = int(query_expression['experiment'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    foreground_decision_func = parse_query_expression(query_expression['foreground'])
    background_decision_func = parse_query_expression(query_expression['background'])
    
    background = [ms for ms in exp.measurements if background_decision_func(ms)]
    foreground = [ms for ms in background if foreground_decision_func(ms)]
    
    enrichment = calculate_feature_enrichment(foreground, background)
    
    motif_tests, motif_results = motif.calculate_motif_enrichment(foreground, background)
    
    background_seqlogo = protein_utils.create_sequence_profile(background)
    background_protein_cnt = len(set( [ms.protein.id for ms in background] ))
    background_peptide_cnt = len( background )
    
    foreground_seqlogo = protein_utils.create_sequence_profile(foreground)
    foreground_protein_cnt = len(set( [ms.protein.id for ms in foreground] ))
    foreground_peptide_cnt = len( foreground )
    
    measurement_data = format_measurement_data(foreground)
    peptide_data = format_peptide_data(foreground)
    
    return {'experiment': exp_id,
            'name': query_expression['name'],
            'background': {
                           'query': query_expression['background'],
                           'proteins': background_protein_cnt,
                           'peptides': background_peptide_cnt,
                           'seqlogo':background_seqlogo
                           },
            'foreground': {
                           'query': query_expression['foreground'],
                           'proteins': foreground_protein_cnt,
                           'peptides': foreground_peptide_cnt,
                           'seqlogo': foreground_seqlogo
                           },
            'measurements': measurement_data,
            'peptides': peptide_data,
            'enrichment': enrichment,
            'motif': {
                      'tests': motif_tests,
                      'results': motif_results
                      }
        }
    

@view_config(route_name='save_subset', renderer='json', permission='private')
def save_subset(request):
    return {}

@view_config(route_name='fetch_subset', renderer='json', permission='private')
def fetch_subset(request):
    return {}