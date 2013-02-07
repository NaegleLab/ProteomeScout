import sys, os, time, csv
from ptmscout.config import strings, settings
from ptmscout.utils import uploadutils


def insert_errors(errors, rows):
    for row in rows:
        row.insert(0, [])
    
    for error in errors:
        rows[error.line-1][0].append(error.message)
        
    for row in rows:
        row[0] = ", ".join(row[0])


def annotate_experiment_file(exp, show_errors):
    header, rows = uploadutils.load_header_and_data_rows(exp.dataset, sys.maxint)
    
    if show_errors:
        header.insert(0, strings.experiment_upload_error_reasons_column_title)
        insert_errors(exp.errors, rows)
        
        rows = [r for r in rows if r[0] != ""]

    return header, rows

def format_experiment_data(exp):
    has_runs = False
    data_headers = set()
    rows = []

    for ms in exp.measurements:
        runs = {}
        for d in ms.data:
            exp_data = runs.get(d.run, [])
            header = (d.type, d.priority, "%s:%s:%s" % (d.type, d.units, d.label))
            data_headers.add(header)

            exp_data.append(d)

            runs[d.run] = exp_data

        modstr = '; '.join([pep.modification.name for pep in ms.peptides])
        row_template = [ms.query_accession, ms.peptide, modstr]

        if len(runs) == 0:
            rows.append(row_template)
        else:
            has_runs = True

        for r in runs:
            sorted_data = sorted(runs[r], key=lambda d: (d.type, d.priority))
            trow = row_template[:] + [r] + [ str(d.value) for d in sorted_data ]
            rows.append(trow)

    final_headers = ['accession', 'peptide', 'modification']
    if has_runs:
        final_headers += ['run']

    data_headers = sorted(list(data_headers), key=lambda item: (item[0], item[1]))
    data_headers = [ d[2] for d in data_headers ]
    final_headers += data_headers

    return final_headers, rows

def experiment_to_tsv(exp):
    headers, rows = format_experiment_data(exp)

    exp_filename = "experiment_data" + str(time.time())
    exp_path = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, exp_filename)

    with open(exp_path, 'w') as exp_file:
        exp_tsv = csv.writer(exp_file, delimiter='\t')
        exp_tsv.writerow(headers)
        for r in rows:
            exp_tsv.writerow(r)

    return headers, exp_filename