from ptmscout.config import strings
from ptmscout.database import protein, modifications
from pyramid.view import view_config

@view_config(route_name='protein_data', renderer='ptmscout:templates/proteins/protein_data.pt')
def protein_experiment_data_view(request):
    pid = int(request.matchdict['id'])
    prot = protein.getProteinById(pid)
    
    experiment_data = {}
    
    mods = modifications.getMeasuredPeptidesByProtein(pid, request.user)
    
    for mod in mods:
        exp_key = (mod.experiment.id, mod.experiment.name)
        ms_data = experiment_data.get(exp_key, [])
        
        phosphopeps = [p.getName() for p in mod.phosphopeps]
        
        run_data = {}
        for d in mod.data:
            values = run_data.get(d.run, set())
            values.add((d.priority, d.label, d.value, d.type))
            run_data[d.run] = values

        sorted_data = []
        for run in run_data:
            sorted_values = [(label, str(value), type_) for (_, label, value, type_) in sorted(run_data[run], key=lambda item: item[0])]
            data_dict = {'run':run, 'values':sorted_values, 'phosphopeps':phosphopeps}
            
            sorted_data.append(data_dict)
            
        sorted_data = sorted(sorted_data, key=lambda item: item['run'])
        
        
        ms_data.extend(sorted_data)
        experiment_data[exp_key] = ms_data
    
    output_data = [ {'id': eid, 'title':name, 'data':experiment_data[(eid,name)]} for (eid, name) in experiment_data if len(experiment_data[(eid,name)]) > 0]
    
    return {'pageTitle': strings.protein_data_page_title,
            'protein': prot,
            'experiment_data':output_data}
    

