from pyramid.view import view_config
from ptmscout.database import experiment, modifications
from ptmscout.config import strings
import base64
import json

def format_pfam_domains(measurements):
    
    domain_map = {}
    for m in measurements:
        prot = m.protein
        
        for d in prot.domains:
            mapset = domain_map.get(d.label, set())
            mapset.add(prot)
            domain_map[d.label] = mapset
            
        if len(prot.domains) == 0:
            mapset = domain_map.get('None', set())
            mapset.add(prot)
            domain_map['None'] = mapset
    
    for domain in domain_map:
        domain_map[domain] = len(domain_map[domain])
    
    domain_list = sorted(domain_map.items(), key=lambda item: -item[1])
    jsondata = base64.b64encode(json.dumps(domain_list))
    
    return {'table':domain_list, 'json':jsondata}

def format_pfam_sites(measurements):
    
    site_map = {}
    
    for m in measurements:
        for p in m.peptides:
            pep = p.peptide
            
            if pep.protein_domain != None:
                pfam_site = pep.protein_domain.label
            else:
                pfam_site = 'None'
            
            mapset = site_map.get(pfam_site, set())
            mapset.add(m)
            site_map[pfam_site] = mapset

    for site in site_map:
        site_map[site] = len(site_map[site])
        
    site_list = sorted(site_map.items(), key=lambda item: -item[1])
    jsondata = base64.b64encode(json.dumps(site_list))
    
    return {'table':site_list, 'json':jsondata}

@view_config(route_name='experiment_pfam', renderer="ptmscout:templates/experiments/experiment_pfam.pt")
def show_pfam_view(request):
    expid = int(request.matchdict['id'])
    exp = experiment.getExperimentById(expid, request.user)
    measurements = modifications.getMeasuredPeptidesByExperiment(expid, request.user)
    
    user_owner = request.user != None and request.user.experimentOwner(exp)
    
    formatted_sites = format_pfam_sites(measurements)
    formatted_domains = format_pfam_domains(measurements)
    
    return {'pageTitle': strings.experiment_pfam_page_title % (exp.name),
            'experiment':exp,
            'sites':formatted_sites,
            'domains':formatted_domains,
            'user_owner': user_owner}