from pyramid.view import view_config
import json
import base64
from ptmscout.database import experiment, annotations
from ptmscout.utils import protein_utils, motif, fisher
from ptmscout.config import strings
from collections import defaultdict

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



def get_valid_subset_labels(experiment_id, user):
    subset_labels = []
    
    for subset in annotations.getSubsetsForUser(experiment_id, user):
        subset_labels.append(subset.name)
    
    return sorted(subset_labels)

def get_user_annotations(exp_id, user, metadata_fields, quantitative_fields, clustering_labels):
    for annotation_set in annotations.getUserAnnotations(exp_id, user):
        values = set()
        for annotation in annotation_set.annotations:
            if annotation.value:
                values.add(annotation.value)
        values = sorted( list(values) ) 
        
        if annotation_set.type == 'cluster':
            clustering_labels[annotation_set.name] = values
        if annotation_set.type == 'nominative':
            metadata_fields[annotation_set.name] = values
        if annotation_set.type == 'numeric':
            quantitative_fields.add(annotation_set.name)

def format_explorer_view(experiment_id, measurements, user):
    quantitative_fields = set()
    
    accessions = []
    metadata_fields = {}
    scansite_fields = {}
    
    cluster_labels = {}
    subset_labels = get_valid_subset_labels(experiment_id, user)
    
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
    
    get_user_annotations(experiment_id, user, metadata_fields, quantitative_fields, cluster_labels)
    
    field_data = {'experiment_id': experiment_id,
                  'quantitative_fields': sorted( list(quantitative_fields) ),
                  'accessions': accessions,
                  'scansite_fields': scansite_fields,
                  'scansite_keys': sorted( scansite_fields.keys() ),
                  'metadata_fields': metadata_fields,
                  'metadata_keys': sorted( metadata_fields.keys() ),
                  'clustering_sets': sorted( list(cluster_labels.keys()) ),
                  'clustering_labels': cluster_labels,
                  'subset_labels': subset_labels}
    
    return {'field_data': base64.b64encode( json.dumps(field_data) )}

class QueryError(Exception):
    pass


def calculate_hypergeometric_enrichment(foreground, background, decision_func):
    foreground_ids = set()
    a = 0
    b = 0
    c = 0
    d = 0
    
    for item in foreground:
        foreground_ids.add(item.id)
        
        if(decision_func(item)):
            a+=1
        else:
            b+=1
    
    for item in background:
        if item.id in foreground_ids:
            continue
        
        if(decision_func(item)):
            c+=1
        else:
            d+=1
        
    pvalue = fisher.FishersExactTest([[a,b],[c,d]]).two_tail_p()
    return pvalue 


def get_proteins(measurements):
    ids = set()
    prots = []
    
    for ms in measurements:
        if ms.protein_id in ids:
            continue
        
        prots.append(ms.protein)
        ids.add(ms.protein_id)
        
    return prots


def calculate_GO_term_enrichment(foreground, background, required_occurences):
    go_term_map = {}
    represented_go_terms = defaultdict(lambda: defaultdict(lambda: 0))
    aspect_map = {'C':'GO-Cellular Component', 'P':"GO-Biological Process", 'F': "GO-Molecular Function"}
    enrichment = []
    
    for prot in foreground:
        for entry in prot.GO_terms:
            t = entry.GO_term
            represented_go_terms[t.aspect][t.id] += 1
            go_term_map[t.id] = t

    def create_has_go_id(go_id):
        def has_go_id(prot):
            return reduce(bool.__or__, [ entry.GO_term.id == go_id for entry in prot.GO_terms ], False)
        return has_go_id
    
    for aspect in represented_go_terms:
        for go_id in represented_go_terms[aspect]:
            if represented_go_terms[aspect][go_id] >= required_occurences:
                p_value = calculate_hypergeometric_enrichment(foreground, background, create_has_go_id(go_id))
                t = go_term_map[go_id]
                enrichment.append(( aspect_map[t.aspect], t.fullName(), p_value ))
            
    return enrichment

def calculate_PfamDomain_enrichment(foreground, background, required_occurences, domain_cutoff):
    enrichment = []
    valid_labels = defaultdict(lambda: 0)
    
    for prot in foreground:
        for dom in prot.domains:
            valid_labels[dom.label] += 1
            
    def create_has_domain(domain_label):
        def has_domain(prot):
            return reduce(bool.__or__, [ dom.label == domain_label and dom.p_value <= domain_cutoff for dom in prot.domains ], False)
        return has_domain
    
    for label in valid_labels:
        if valid_labels[dom.label] >= required_occurences:
            p_value = calculate_hypergeometric_enrichment(foreground, background, create_has_domain(label))
            enrichment.append( ( 'Pfam-Domain', label, p_value ) )
    
    return enrichment

def calculate_PfamSite_enrichment(foreground, background, required_occurences, domain_cutoff):
    enrichment = []
    valid_labels = defaultdict(lambda: 0)
    
    for ms in foreground:
        for modpep in ms.peptides:
            if modpep.peptide.protein_domain:
                valid_labels[modpep.peptide.protein_domain.label] += 1
                
    def create_has_domain(domain_label):
        def has_domain(ms):
            return reduce(bool.__or__, [ modpep.peptide.protein_domain != None and modpep.peptide.protein_domain.label == domain_label and modpep.peptide.protein_domain.p_value <= domain_cutoff for modpep in ms.peptides ], False)
        return has_domain
    
    for label in valid_labels:
        if valid_labels[label] >= required_occurences:
            p_value = calculate_hypergeometric_enrichment(foreground, background, create_has_domain(label))
            enrichment.append( ( 'Pfam-Site', label, p_value ))
    
    return enrichment

def calculate_Scansite_enrichment(foreground, background, required_occurences, scansite_cutoff):
    enrichment = []
    valid_labels = defaultdict(lambda: defaultdict(lambda: 0))
    source_map = {'scansite_bind': "Scansite-Bind", 'scansite_kinase': "Scansite-Kinase"}
    for ms in foreground:
        for modpep in ms.peptides:
            for prediction in modpep.peptide.predictions:
                valid_labels[prediction.source][prediction.value] += 1
                
    def create_has_scansite(source, label):
        def has_scansite(ms):
            return reduce(bool.__or__, [ prediction.percentile <= scansite_cutoff and prediction.source == source and prediction.value == label for modpep in ms.peptides for prediction in modpep.peptide.predictions ], False)
        return has_scansite        
        
    for source in valid_labels:
        for label in valid_labels[source]:
            if valid_labels[source][label] >= required_occurences:
                p_value = calculate_hypergeometric_enrichment(foreground, background, create_has_scansite(source, label))
                enrichment.append( (source_map[source], label, p_value) )

    return enrichment


def calculate_Region_enrichment(foreground, background, required_occurences):
    enrichment = []
    valid_labels = defaultdict(lambda: 0)
    for ms in foreground:
        for region in ms.protein.regions:
            valid_labels[region.label] += 1
            
    def create_has_region(label):
        def has_region(ms):
            return reduce(bool.__or__, [ region.hasSite( modpep.peptide.site_pos ) for region in ms.protein.regions for modpep in ms.peptides ], False)
        return has_region
    
    for label in valid_labels:
        if valid_labels[label] >= required_occurences:
            p_value = calculate_hypergeometric_enrichment(foreground, background, create_has_region(label))
            enrichment.append( ( 'Region', label, p_value ) )
    
    return enrichment


def calculate_annotation_enrichment(name, foreground, background, required_occurences):
    enrichment = []
    
    valid_labels = defaultdict(lambda: 0)
    for ms in foreground:
        if name in ms.annotations and ms.annotations[name] != None:
            valid_labels[ms.annotations[name]] += 1
    
    def create_has_annotation(label):
        def has_annotation(ms):
            return name in ms.annotations and ms.annotations[name] == label
        return has_annotation
    
    for label in valid_labels:
        if valid_labels[label] >= required_occurences:
            p_value = calculate_hypergeometric_enrichment(foreground, background, create_has_annotation(label))
            enrichment.append( ( "Annotation: %s" % (name), label, p_value ) )
    
    return enrichment


def calculate_feature_enrichment(foreground, background, annotation_types, required_occurences=1, scansite_cutoff=5, domain_cutoff=1):
    enrichment_table = []

    p_foreground = get_proteins(foreground)
    p_background = get_proteins(background)
    
    enrichment_table += calculate_GO_term_enrichment(p_foreground, p_background, required_occurences)
    enrichment_table += calculate_PfamDomain_enrichment(p_foreground, p_background, required_occurences, domain_cutoff)
    
    enrichment_table += calculate_PfamSite_enrichment(foreground, background, required_occurences, domain_cutoff)
    enrichment_table += calculate_Scansite_enrichment(foreground, background, required_occurences, scansite_cutoff)
    enrichment_table += calculate_Region_enrichment(foreground, background, required_occurences)
    
    for aset in annotation_types:
        if annotation_types[aset] == 'nominative':
            enrichment_table += calculate_annotation_enrichment(aset, foreground, background, required_occurences)
    
    return sorted(enrichment_table, key=lambda item: (item[2], item[0], item[1]) )

def create_protein_filter(filter_value):
    def protein_filter(ms):
        return ms.protein.hasAccession(filter_value) or \
                ms.protein.name == filter_value or \
                ms.protein.acc_gene == filter_value or \
                ms.protein.locus == filter_value
    return protein_filter


def create_metadata_filter(field, op, value, annotation_types):
    value = value.strip()
    
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
        return [ None != modpep.peptide.protein_domain and modpep.peptide.protein_domain.label.lower() == value.lower() for modpep in ms.peptides ]
    
    @op_apply_wrapper
    def region_filter(ms):
        return [ region.label.lower() == value.lower() and region.hasSite(modpep.peptide.site_pos) for region in ms.protein.regions for modpep in ms.peptides ]
                
    def create_custom_annotation_filter(field):
        def custom_filter(ms):
            if (field not in ms.annotation_types) or (None == ms.annotations[field]):
                return []
            if ms.annotations[field].lower() == value.lower():
                return [ True ]
            
        return custom_filter
    
    if field == 'Modification Type':
        return modtype_filter
    elif field == 'Modified Residue':
        return modsite_filter
    elif field.find('GO-') == 0:
        return GO_filter
    elif field == 'Pfam-Domain':
        return pfam_domain_filter
    elif field == 'Pfam-Site':
        return pfam_site_filter
    elif field == 'Region':
        return region_filter
    elif field in annotation_types and annotation_types[field] in set(['cluster','nominative']):
        return create_custom_annotation_filter(field)

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
        return [prediction.source == 'scansite_kinase' and prediction.value == value and prediction.percentile < stringency 
                 for modpep in ms.peptides for prediction in modpep.peptide.predictions]
    
    @op_apply_wrapper
    def bind_filter(ms):
        return [prediction.source == 'scansite_bind' and prediction.value == value and prediction.percentile < stringency 
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

def parse_value(value, annotation_types):
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
                
            for aset in ms.annotation_types:
                if aset == value and ms.annotation_types[aset] == 'numeric':
                    avalue = ms.annotations[aset]
                    if avalue == None:
                        raise MissingMeasurementError()
                    return float(avalue)
            
            if not (value in annotation_types and annotation_types[value] == 'numeric'):
                raise QueryError("No such data label: '%s'" % (value))
            
            raise MissingMeasurementError()

        return get_data

def parse_arthmetic_expression(expression, annotation_types):
    if len(expression) == 1:
        return parse_value(expression[0], annotation_types)
    elif len(expression) == 3:
        v1 = parse_value(expression[0], annotation_types)
        op = expression[1]
        v2 = parse_value(expression[2], annotation_types)
        
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
        

def create_quantitative_filter(expression, annotation_types):
    comparison_operators = set([ 'eq', 'neq', 'gt', 'lt', 'geq', 'leq' ])
    
    comp_index = find_in(expression, comparison_operators)
    if len(comp_index) == 0:
        raise QueryError("No comparison found in quantitative expression: %s" % (readable_expression(expression)))
    comp_index = comp_index[0]
    
    LHS = expression[:comp_index]
    RHS = expression[comp_index+1:]
    
    LHS_func = parse_arthmetic_expression(LHS, annotation_types)
    RHS_func = parse_arthmetic_expression(RHS, annotation_types)
    
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


def create_subset_filter(op, subset_label, annotation_types, experiment_id=None, user=None):
    subset = annotations.getSubsetByName(experiment_id, subset_label, user)
    if subset == None:
        raise QueryError("No such subset '%s'" % (subset_label))
    
    return parse_query_expression(subset.foreground_query, experiment_id, user, annotation_types)
  

def parse_condition(condition, experiment_id, user, annotation_types):
    if condition[0] == 'protein':
        return create_protein_filter(condition[1])
    if condition[0] == 'metadata':
        return create_metadata_filter(*tuple(condition[1:] + [ annotation_types ]))
    if condition[0] == 'scansite':
        return create_scansite_filter(*tuple(condition[1:]))
    if condition[0] == 'sequence':
        return create_sequence_filter(condition[1])
    if condition[0] == 'quantitative':
        return create_quantitative_filter(condition[1:], annotation_types)
    if condition[0] == 'subset':
        return create_subset_filter(*tuple(condition[1:] + [ annotation_types ]), user=user, experiment_id=experiment_id)
    if condition[0] == 'cluster':
        return create_metadata_filter(*tuple(condition[1:] + [ annotation_types ]))

    raise QueryError("No filter type found for '%s'" % (condition[0]))

def parse_dnf_clause(conditions, experiment_id, user, annotation_types):
    def conjunction(operands):
        def apply_conjunction(ms):
            for op in operands:
                if not op(ms):
                    return False
            return True
        
        return apply_conjunction
    
    return conjunction( [ parse_condition(cond, experiment_id, user, annotation_types) for cond in conditions] )

def parse_query_expression(expression, experiment_id, user, annotation_types):
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
        
    return disjunction( [ parse_dnf_clause(op, experiment_id, user, annotation_types) for op in dnf_clauses ] )

def get_peptide_annotations(peptide, protein):
    domains = [ d for d in protein.domains if d.start <= peptide.site_pos and peptide.site_pos <= d.stop ]
    regions = [ r for r in protein.regions if r.start <= peptide.site_pos and peptide.site_pos <= r.stop ]
    mutations = [ m for m in protein.mutations if m.location == peptide.site_pos ]

    annotations = {'kinase':False, 'loop':None, 'mutant':False}

    if len(domains) > 0:
        annotations['kinase'] = domains[0].label.find('Pkinase') >= 0
    if len(regions) > 0:
        if regions[0].label == 'Kinase Activation Loop':
            annotations['loop'] = 'K'
        elif regions[0].label == 'Possible Kinase Activation Loop':
            annotations['loop'] = '?'
    if len(mutations) > 0:
        annotations['mutant'] = True

    return annotations

def format_peptide_data(request, experiment_id, measurements):
    formatted_peptide_data = []
    for ms in measurements:
        gene_name = ms.protein.getGeneName()
        protein_id = ms.protein.id
        protein_link = request.route_url('protein_main', id=protein_id) + "?experiment_id=%d" % (experiment_id)
        
        for modpep in ms.peptides:
            pep_aligned = modpep.peptide.pep_aligned
            site_name = modpep.peptide.getName()
            mod_name = modpep.modification.name
            
            annotations = get_peptide_annotations(modpep.peptide, ms.protein)
            formatted_peptide_data.append( (ms.id, ms.query_accession, gene_name, protein_link, pep_aligned, site_name, mod_name, annotations) )
    
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

def compute_annotations(exp, user, measurements):
    annotation_sets = annotations.getUserAnnotations(exp.id, user)
    annotation_types = {}
    ms_map = {}
    for ms in measurements:
        ms_map[ms.id] = ms
        ms.annotations = {}
        ms.annotation_types = {}
    
    for aset in annotation_sets:
        annotation_types[aset.name] = aset.type
        for annotation in aset.annotations:
            ms = ms_map[annotation.MS_id]
            ms.annotations[aset.name] = annotation.value
            ms.annotation_types[aset.name] = aset.type
            
    return annotation_types

def compute_subset_enrichment(request, exp, user, subset_name, foreground_exp, background_exp):
    measurements = exp.measurements
    annotation_types = compute_annotations(exp, user, measurements)
    
    foreground_decision_func = parse_query_expression(foreground_exp, exp.id, user, annotation_types)
    background_decision_func = parse_query_expression(background_exp, exp.id, user, annotation_types)
    
    background = [ms for ms in measurements if background_decision_func(ms)]
    foreground = [ms for ms in background if foreground_decision_func(ms)]
    
    if len(foreground) == 0:
        return {'status':"error",
                'message':strings.error_message_subset_empty}
    
    enrichment = calculate_feature_enrichment(foreground, background, annotation_types)
    
    motif_tests, motif_results = motif.calculate_motif_enrichment(foreground, background)
    
    background_seqlogo = protein_utils.create_sequence_profile(background)
    background_site_cnt = len(set([(modpep.peptide_id, modpep.modification_id) for ms in background for modpep in ms.peptides]))
    background_protein_cnt = len(set( [ms.protein.id for ms in background] ))
    background_peptide_cnt = len( background )
    
    foreground_seqlogo = protein_utils.create_sequence_profile(foreground)
    foreground_site_cnt = len(set([(modpep.peptide_id, modpep.modification_id) for ms in foreground for modpep in ms.peptides]))
    foreground_protein_cnt = len(set( [ms.protein.id for ms in foreground] ))
    foreground_peptide_cnt = len( foreground )
    
    measurement_data = format_measurement_data(foreground)
    peptide_data = format_peptide_data(request, exp.id, foreground)
    
    return {'status':'success',
            'experiment': exp.id,
            'name': subset_name,
            'background': {
                           'query': background_exp,
                           'proteins': background_protein_cnt,
                           'peptides': background_peptide_cnt,
                           'sites': background_site_cnt,
                           'seqlogo':background_seqlogo
                           },
            'foreground': {
                           'query': foreground_exp,
                           'proteins': foreground_protein_cnt,
                           'peptides': foreground_peptide_cnt,
                           'sites': foreground_site_cnt,
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

@view_config(route_name='compute_subset', request_method='POST', renderer='json', permission='private')
def perform_subset_enrichment(request):
    query_expression = json.loads(request.body, encoding=request.charset)
    exp_id = int(query_expression['experiment'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    subset_name = query_expression['name']
    foreground_exp = query_expression['foreground']
    background_exp = query_expression['background']
    
    return compute_subset_enrichment(request, exp, request.user, subset_name, foreground_exp, background_exp)
    

@view_config(route_name='save_subset', request_method='POST', renderer='json', permission='private')
def save_subset(request):
    query_expression = json.loads(request.body, encoding=request.charset)
    exp_id = int(query_expression['experiment'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    subset_name = query_expression['name']
    foreground_exp = query_expression['foreground']
    background_exp = query_expression['background']
    
    annotation_types = compute_annotations(exp, request.user, exp.measurements)
    parse_query_expression(foreground_exp, exp_id, request.user, annotation_types)
    parse_query_expression(background_exp, exp_id, request.user, annotation_types)
    
    if(None != annotations.getSubsetByName(exp_id, subset_name, request.user)):
        return {'status':"error",
                'message':strings.error_message_subset_name_exists}
    
    subset = annotations.Subset()
    subset.name = subset_name
    subset.experiment_id = exp_id
    subset.owner_id = request.user.id
    subset.foreground_query = foreground_exp
    subset.background_query = background_exp
    subset.save()
    
    return {'id': subset.id,
            'name': subset.name,
            'status':'success'}

@view_config(route_name='fetch_subset', request_method='POST', renderer='json', permission='private')
def fetch_subset(request):
    query_expression = json.loads(request.body, encoding=request.charset)
    exp_id = int(query_expression['experiment'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    subset = None
    if 'id' in query_expression:
        subset = annotations.getSubsetById(int(query_expression['id']), exp_id)
    elif 'name' in query_expression:
        subset = annotations.getSubsetByName(exp_id, query_expression['name'], request.user)
    
    if None == subset:
        return {'status':"error",
                'message':strings.error_message_subset_name_does_not_exist}
    
    
    result = compute_subset_enrichment(request, exp, request.user, subset.name, subset.foreground_query, subset.background_query)
    result['id'] = subset.id
    return result
