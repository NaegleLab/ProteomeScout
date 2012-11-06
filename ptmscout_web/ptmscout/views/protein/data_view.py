from ptmscout.config import strings
from ptmscout.database import protein, modifications
from pyramid.view import view_config
from ptmscout.utils import webutils

def format_protein_data(mods):
    experiment_data = {}
    
    for mod in mods:
        exp_key = (mod.experiment.id, mod.experiment.name)
        ms_data = experiment_data.get(exp_key, [])
        
        phosphopeps = [p.getName() for p in mod.phosphopeps]
        
        run_data = {}
        for d in mod.data:
            values = run_data.get(d.run, set())
            values.add((d.priority, d.units, d.label, d.value, d.type))
            run_data[d.run] = values

        sorted_data = []
        for run in run_data:
            data_units = [item[1] for item in run_data[run] if item[4] == 'data']
            units = data_units[0]
            
            sorted_values = [(label, str(value), type_) for (_, _, label, value, type_) in sorted(run_data[run], key=lambda item: item[0])]
            data_dict = {'run':run, 'units': units, 'values':sorted_values, 'phosphopeps':phosphopeps}
            
            sorted_data.append(data_dict)
            
        sorted_data = sorted(sorted_data, key=lambda item: item['run'])
        
        
        ms_data.extend(sorted_data)
        experiment_data[exp_key] = ms_data
    
    return [ {'id': eid, 'title':name, 'data':experiment_data[(eid,name)]} for (eid, name) in experiment_data if len(experiment_data[(eid,name)]) > 0]

@view_config(route_name='protein_data', renderer='ptmscout:templates/proteins/protein_data.pt')
def protein_experiment_data_view(request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    experiment_id = webutils.get(request, 'experiment_id', None)
    
    if experiment_id == None:
        mods = modifications.getMeasuredPeptidesByProtein(pid, request.user)
    else:
        mods = modifications.getMeasuredPeptidesByExperiment(int(experiment_id), request.user, [pid])
    
    output_data = format_protein_data(mods)
        
    return {'pageTitle': strings.protein_data_page_title,
            'protein': prot,
            'experiment_data':output_data}
    

