import sys, os, csv, time, pickle
from ptmscout.config import settings

def format_size( size ):
    postfix = ['B','KB','MB','GB','TB']

    i = 0
    while size > 1024:
        size /= float(1024)
        i+=1

    if type(size) == long:
        return "%d %s" % ( size, postfix[i] )
    else:
        return "%.1f %s" % ( size, postfix[i] )

if __name__=='__main__':
    summary_struct = {}
    for filename in sys.argv[1:]:
        if filename.find('.tsv') == -1:
            continue
        fpath = os.path.join(settings.ptmscout_path, settings.export_file_path, filename)

        i = 0
        j = 0
        with open(fpath, 'r') as f:
            dr = csv.DictReader(f, dialect='excel-tab')

            for row in dr:
                i+=1
                mods = row['modifications'].strip().split(';')
                if len(mods) == 1 and mods[0].strip() == '':
                    mods.pop()
                j += len(mods)

        mod_time = time.ctime(os.path.getmtime(fpath))
        size = format_size( os.path.getsize(fpath) )
        summary_struct[filename] = {'proteins':i,'modifications':j,'date':mod_time,'size':size}

    pickle.dump(summary_struct, sys.stdout)
