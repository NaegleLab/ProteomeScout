from pyramid.view import view_config
from ptmscout.utils import webutils
from ptmscout.database import experiment, modifications
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.utils import downloadutils
from ptmscout.config import settings

@view_config(route_name='experiment_download', renderer='tsv', permission='private')
def download_experiment(request):
    get_errors = bool( webutils.get(request, 'errors', False) )
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    if exp not in request.user.myExperiments():
        raise HTTPForbidden()

    header, rows = downloadutils.annotate_experiment_with_errors(exp, get_errors)
    
    exp_filename = 'experiment.%d.tsv' % (exp_id)
    if get_errors:
        exp_filename = 'experiment.%d.errors.tsv' % (exp_id)
   
    request.response.content_type = 'text/tab-separated-values'
    request.response.content_disposition = 'attachment; filename="%s"' % (exp_filename)
    return { 'header': header, 'data': rows }

def annotate_experiment(user, exp, header, rows):
    header += ['nearby_modifications', 'nearby_mutations', 'site_domains', 'site_regions']
    protein_mods = {}
    
    for ms in exp.measurements:
        if ms.protein_id not in protein_mods:
            protein_mods[ms.protein_id] = modifications.getMeasuredPeptidesByProtein(ms.protein_id, user)
    
    ms_map = {}    
    for ms in exp.measurements:
        ms_map[ms.id] = ms
        
    for row in rows:
        ms = ms_map[row[0]]
        prot = ms.protein
        
        min_range = ms.peptides[0].peptide.site_pos - 7
        max_range = ms.peptides[-1].peptide.site_pos + 7
        
        nearby_modifications = set()
        for ms2 in protein_mods[ms.protein_id]:
            for modpep in ms2.peptides:
                site_type = modpep.peptide.site_type
                site_pos = modpep.peptide.site_pos
                mod_name = modpep.modification.name
                
                if min_range <= site_pos and site_pos <= max_range:
                    nearby_modifications.add((site_pos, site_type, mod_name))
                    
        nearby_modifications = [ "%s%d: %s" % (site_type, site_pos, mod_name) for site_pos, site_type, mod_name in sorted(list(nearby_modifications)) ]
        nearby_mutations = [ str(mutation) for mutation in sorted(prot.mutations, key=lambda item: item.location) if min_range < mutation.location and mutation.location < max_range ]
        
        site_domains = list(set([ modpep.peptide.protein_domain.label for modpep in ms.peptides if modpep.peptide.protein_domain ]))
        site_regions = list(set([ region.label for modpep in ms.peptides for region in prot.regions if region.hasSite(modpep.peptide.site_pos) ]))

        sep = settings.mod_separator_character + ' '

        row.append(sep.join(nearby_modifications))
        row.append(sep.join(nearby_mutations))
        row.append(sep.join(site_domains))
        row.append(sep.join(site_regions))
        
def get_experiment_header(exp):
    header = ['MS_id', 'query_accession', 'peptide', 'mod_sites', 'aligned_peptides', 'modification_types']
    
    data_labels = set()
    for ms in exp.measurements:
        for d in ms.data:
            data_labels.add((d.run,d.type,d.units,d.label))
    
    def float_last_term(r,dt,u,l):
        try:
            l = float(l)
        except:
            pass
        
        return (r,dt,u,l)
    
    data_labels = [ "%s:%s:%s:%s" % (r,dt,u,str(l)) for r,dt,u,l in sorted(list(data_labels), key=lambda item: float_last_term(*item)) ]
    header += data_labels
    
    return header, data_labels

def get_experiment_data(exp, data_labels):
    rows = []
    for ms in exp.measurements:
        mod_sites = '; '.join([modpep.peptide.getName() for modpep in ms.peptides])
        aligned_peptides = '; '.join([modpep.peptide.pep_aligned for modpep in ms.peptides])
        modification_types = '; '.join([modpep.modification.name for modpep in ms.peptides])
        
        row = [ms.id, ms.query_accession, ms.peptide, mod_sites, aligned_peptides, modification_types]
        
        ms_data = {}
        for d in ms.data:
            formatted_label = "%s:%s:%s:%s" % (d.run, d.type, d.units, d.label)
            ms_data[formatted_label] = d.value
            
        for dl in data_labels:
            row.append( ms_data[dl] )
            
        rows.append(row)
            
    return rows

@view_config(route_name='experiment_export', renderer='tsv', permission='private')
def export_experiment(request):
    annotate = webutils.get(request, 'annotate', 'no') == 'yes'
    exp_id = int(request.matchdict['id'])
    exp = experiment.getExperimentById(exp_id, request.user)
    
    header, data_labels = get_experiment_header(exp)
    rows = get_experiment_data(exp, data_labels)
    
    exp_filename = 'experiment.%d.export.tsv' % (exp_id)
    if annotate:
        exp_filename = 'experiment.%d.annotated.tsv' % (exp_id)
        annotate_experiment(request.user, exp, header, rows)    
    
    request.response.content_type = 'text/tab-separated-values'
    request.response.content_disposition = 'attachment; filename="%s"' % (exp_filename)
    return { 'header': header, 'data': rows }