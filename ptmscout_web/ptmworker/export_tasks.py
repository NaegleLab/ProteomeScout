from ptmworker.helpers import upload_helpers
from ptmscout.config import settings, strings
from ptmscout.database import experiment, modifications, user, modifications
import celery
from ptmworker import notify_tasks
import csv, os

def annotate_experiment(user, exp, header, rows):
    header += [ 'scansite_bind', 'scansite_kinase', 'nearby_modifications', 'nearby_mutations', 'site_domains', 'site_regions',\
                    'protein_domains', 'protein_GO_BP', 'protein_GO_CC', 'protein_GO_MF' ]
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
        
        scansite_kinase = []
        scansite_bind = []
        for modpep in ms.peptides:
            for pp in modpep.peptide.predictions:
                if pp.source=='scansite_kinase':
                    scansite_kinase.append( "%s (%.2f)" % ( pp.value, pp.percentile ))
                if pp.source=='scansite_bind':
                    scansite_bind.append( "%s (%.2f)" % ( pp.value, pp.percentile ))

        protein_domains = []
        protein_GO = { 'P':set(), 'F':set(), 'C':set() }

        for d in prot.domains:
            protein_domains.append( "%s (%d-%d)" % ( d.label, d.start, d.stop ) )
        for ge in prot.GO_terms:
            term = ge.GO_term
            protein_GO[term.aspect].add(term.GO)

        row.append(sep.join(scansite_bind))
        row.append(sep.join(scansite_kinase))

        row.append(sep.join(nearby_modifications))
        row.append(sep.join(nearby_mutations))
        row.append(sep.join(site_domains))
        row.append(sep.join(site_regions))

        row.append(sep.join(protein_domains))
        row.append(sep.join(list(protein_GO['P'])))
        row.append(sep.join(list(protein_GO['C'])))
        row.append(sep.join(list(protein_GO['F'])))
        
def get_experiment_header(exp):
    header = ['MS_id', 'query_accession', 'gene', 'locus', 'protein_name', 'species', 'peptide', 'mod_sites', 'gene_site', 'aligned_peptides', 'modification_types']
    
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
        
        gene_sites = [ms.protein.getGeneName()] + [modpep.peptide.getName() for modpep in ms.peptides]
        row = [ms.id, ms.query_accession, ms.protein.acc_gene, ms.protein.locus, ms.protein.name, ms.protein.species.name, ms.peptide, mod_sites, '_'.join(gene_sites), aligned_peptides, modification_types]
        
        ms_data = {}
        for d in ms.data:
            formatted_label = "%s:%s:%s:%s" % (d.run, d.type, d.units, d.label)
            ms_data[formatted_label] = d.value
            
        for dl in data_labels:
            row.append( ms_data[dl] )
            
        rows.append(row)
            
    return rows


@celery.task
@upload_helpers.notify_job_failed
@upload_helpers.logged_task
def run_experiment_export_job(annotate, export_id, exp_id, user_id, job_id):
    notify_tasks.set_job_status.apply_async((job_id, 'started'))
    notify_tasks.set_job_stage.apply_async((job_id, 'exporting', 0))

    exp_filename = 'experiment.%d.%d.%d.tsv' % (exp_id, user_id, export_id)
    exp_path = os.path.join(settings.ptmscout_path, settings.annotation_export_file_path, exp_filename)

    usr = user.getUserById(user_id)
    exp = experiment.getExperimentById(exp_id, usr)

    header, data_labels = get_experiment_header(exp)
    rows = get_experiment_data(exp, data_labels)

    if annotate:
        annotate_experiment(usr, exp, header, rows)

    with open(exp_path, 'w') as export_file:
        cw = csv.writer(export_file, dialect='excel-tab')

        cw.writerow(header)
        for row in rows:
            cw.writerow(row)

    notify_tasks.finalize_experiment_export_job.apply_async((job_id,))
