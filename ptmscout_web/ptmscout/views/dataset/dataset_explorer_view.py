from pyramid.view import view_config
import json
import base64

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

def get_protein_metadata(measurements, metadata_fields):
    metadata_fields.update({
                       'Gene':set(),
                       'Protein Accession':set(),
                       'Protein Name':set(),
                       'Species':set(),
                       })
    
    for ms in measurements:
        metadata_fields['Gene'].add(ms.protein.getGeneName())
        metadata_fields['Protein Name'].add(ms.protein.name)
        metadata_fields['Species'].add(ms.protein.species.name)
        
        for acc in ms.protein.accessions:
            metadata_fields['Protein Accession'].add(acc.value)
            
        
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



def get_GO_metadata(measurements, metadata_fields):
    metadata_fields.update({'GO-Cellular Component':set(),
                            'GO-Molecular Function':set(),
                            'GO-Biological Process':set()})

    for ms in measurements:
        for go_term in ms.protein.GO_terms:
            if go_term.GO_term.aspect == 'C':
                metadata_fields['GO-Cellular Component'].add(go_term.GO_term.term)
            if go_term.GO_term.aspect == 'P':
                metadata_fields['GO-Biological Process'].add(go_term.GO_term.term)
            if go_term.GO_term.aspect == 'F':
                metadata_fields['GO-Molecular Function'].add(go_term.GO_term.term)


def get_modification_metadata(measurements, metadata_fields):
    metadata_fields.update({'Modified Residue':set(),
                            'Modification Type':set()})

    for ms in measurements:
        for modpep in ms.peptides:
            metadata_fields['Modified Residue'].add(modpep.peptide.site_type)
            
            metadata_fields['Modification Type'].add(modpep.modification.name)
            for ptm in modpep.modification.getAllParents():
                metadata_fields['Modification Type'].add(ptm.name)


def format_explorer_view(measurements):
    quantitative_fields = set()
    
    metadata_fields = {}
    cluster_labels = {}
    subset_labels = []
    
    get_protein_metadata(measurements, metadata_fields)
    get_pfam_metadata(measurements, metadata_fields)
    get_scansite_metadata(measurements, metadata_fields)
    get_GO_metadata(measurements, metadata_fields)
    get_modification_metadata(measurements, metadata_fields)
    
    for ms in measurements:
        for d in ms.data:
            quantitative_fields.add("%s:%s:%s" % (d.run, d.type, d.label))

    
    for field in metadata_fields.keys():
        metadata_fields[field] = list(metadata_fields[field])
    
    field_data = {'quantitative_fields': sorted( list(quantitative_fields) ),
                  'metadata_fields': metadata_fields,
                  'metadata_keys': sorted( metadata_fields.keys() ),
                  'clustering_sets': sorted( list(cluster_labels.keys()) ),
                  'clustering_labels': cluster_labels,
                  'subset_labels':subset_labels}
    
    return {'field_data': base64.b64encode( json.dumps(field_data) )}



@view_config(route_name='compute_subset', renderer='json', permission='private')
def compute_subset_enrichment(request):
    return {}