from ptmscout.config import settings
import os
import time
from ptmscout.utils import decorators

def parse_motif_result(result_file, log_file):
    
    lf = open(log_file, 'r')
    log_values = {}
    for line in lf:
        k,v = line.strip().split()
        log_values[k] = v
    lf.close()
    
    rf = open(result_file, 'r')
    results = []
    for line in rf:
        _, pep, fgf, bgf, pv, _ = line.strip().split("|")
        
        fg_matches, fg_total = fgf.split('/')
        bg_matches, bg_total = bgf.split('/')
        
        pv = float(pv.strip())
        
        results.append((pep.strip(), (int(fg_matches), int(fg_total)), (int(bg_matches), int(bg_total)), pv))
    
    rf.close()
        
    return int(log_values["NUM_TESTS_FOR_KRISTEN"]), results
        
    
def save_peptides(measurements, filename):
    ff = open(filename, 'w')
    for ms in measurements:
        for mspep in ms.peptides:
            ff.write("%s\n" % (mspep.peptide.pep_aligned))
    ff.close()

def run_motif_result(foreground_file, background_file, result_file):
    os.system('./run_motif_enrichment.sh %s %s %s' % (foreground_file, background_file, result_file))


@decorators.pushdir( os.path.join(settings.ptmscout_path, settings.motif_script_path) )
def calculate_motif_enrichment(foreground, background):
    postfix = str(time.time())
    foreground_file = "foreground%s" % (postfix)
    background_file = "background%s" % (postfix)
    result_file = "result%s" % (postfix)
    result_log = "result%s.log" % (postfix)
    
    save_peptides(foreground, foreground_file)
    save_peptides(background, background_file)
    
    run_motif_result(foreground_file, background_file, result_file)
    test_cnt, results = parse_motif_result(result_file, result_log)
    
    os.remove(foreground_file)
    os.remove(background_file)
    os.remove(result_file)
    os.remove(result_log)
    
    return test_cnt, results