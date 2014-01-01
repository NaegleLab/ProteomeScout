from ptmscout.config import settings
import os
import time
from ptmscout.utils import decorators
import random

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
        
    return int(log_values["NUM_TESTS_FOR_KRISTEN"]), sorted(results, key=lambda item: item[3])
    
    
def save_peptides(measurements, filename):
    ff = open(filename, 'w')
    keys_seen = set()
    for ms in measurements:
        for mspep in ms.peptides:
            pepkey = (mspep.peptide_id, mspep.modification_id)
            if pepkey not in keys_seen:
                ff.write("%s\n" % (mspep.peptide.pep_aligned))
                keys_seen.add(pepkey)
    ff.close()

def run_motif_result(foreground_file, background_file, result_file, k):
    command = '%s %s %s %s --length %d' % (os.path.join(settings.ptmscout_path, 'scripts', 'motif', 'run_motif_enrichment.sh'), foreground_file, background_file, result_file, k)
    os.system(command)

#@decorators.pushdir( os.path.join(settings.ptmscout_path, settings.motif_script_path) )
def calculate_motif_enrichment(foreground, background, k = 5):
    postfix = str(time.time())
    random_postfix = random.randint(0,1000000)
    
    scripts_dir = os.path.join(settings.ptmscout_path, 'scripts', 'motif')
    
    foreground_file = "foreground%s.%d" % (postfix, random_postfix)
    foreground_path = os.path.join(scripts_dir, foreground_file)
    background_file = "background%s.%d" % (postfix, random_postfix)
    background_path = os.path.join(scripts_dir, background_file)
    result_file = os.path.join(scripts_dir, "result%s.%d" % (postfix, random_postfix))
    result_log = os.path.join(scripts_dir, "result%s.%d.log" % (postfix, random_postfix))
    
    save_peptides(foreground, foreground_path)
    save_peptides(background, background_path)
    
    run_motif_result(foreground_file, background_file, result_file, k)
    test_cnt, results = parse_motif_result(result_file, result_log)
    
    os.remove(foreground_path)
    os.remove(background_path)
    os.remove(result_file)
    os.remove(result_log)
    
    return test_cnt, results
