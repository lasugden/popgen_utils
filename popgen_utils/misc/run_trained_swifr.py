import multiprocessing as mp
import subprocess

from datetime import datetime
import os
import os.path as opath
import yaml
try:
    import importlib.resources as ilresources
except ImportError:
    try:
        import importlib_resources as ilresources
    except ImportError:
        raise ImportError('Must install backport of importlib_resources if not using Python >= 3.7')

from popgen_utils import config
from popgen_utils.misc import hashing
from popgen_utils.misc import oscar_scripts

def run_swifr_on_training_data(project_name, model_name, type, pop_of_interest, 
    swifr_trained_path, out_path, pi_vec, data_path=None):
    """
    Iterate over populations, selection strengths, and sweep times and extract combined statistics.

    Args:
        project name (str): name of project (directory)
        model_name (str): e.g. gravel_neutral_400kb
        type (str): 'sweep' or 'neutral'
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        swifr_trained_path (str): path to be passed to swifr_test where trained model lives (relative path, from common
        data directory)
        out_path (str): path for the output (also relative path from common data directory)
        pi_vec (array): pi values for neutral, sweep, linked (in that order)

    """
    pops = ['p1','p2','p3']
    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path, model_name)
        bash_path = opath.join(slim_model_path,'bash')

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)

    txt = ilresources.open_text(oscar_scripts, 'run_trained_swifr.txt').read()
    
    if type == 'sweep':
        for scoeff in params['selection_coefficient']:
            for time in params['sweep_time']:
                pops_reference = [x for x in pops if x != pop_of_interest]
                parameter_model_name = (f'{model_name}_coeff-{scoeff}_'
                                            f'pop-{pop_of_interest}_start-{time}')
                allstats_file = opath.join(slim_model_path, 
                    parameter_model_name+'_%s_%s_allstats.txt' % (pop_of_interest, ''.join(pops_reference)))
                path2trained = opath.join(slim_path, swifr_trained_path)
                outfile = opath.join(slim_path, out_path, 'sweep', parameter_model_name+'_classified')
                if not opath.exists(outfile):
                    formatted_txt = txt.format(**{
                        'path2trained' : path2trained,
                        'allstats_file' : allstats_file,
                        'pi_vals' : ' '.join([str(x) for x in pi_vec]),
                        'outfile' : outfile,
                        'log_file' : opath.join(bash_path,parameter_model_name+'_classify_'+pop_of_interest),
                        })
                    fp = open(opath.join(bash_path, parameter_model_name+'_classify_'+pop_of_interest+'.sh'),'w')
                    fp.write(formatted_txt)
                    fp.close()
                    os.system('sbatch '+opath.join(bash_path, parameter_model_name+'_classify_'+pop_of_interest+'.sh'))
                    #os.system('swifr_test --path2trained '+path2trained+' --file '+allstats_file+
                    #   ' --pi '+' '.join([str(x) for x in pi_vec])+' --outfile '+outfile)

    elif type == 'neutral':
        for sim in params['sims']:
            pops_reference = [x for x in pops if x!= pop_of_interest]
            parameter_model_name = (f'{model_name}_sim-{sim}')
            allstats_file = opath.join(slim_model_path, 
                parameter_model_name+'_%s_%s_allstats.txt' % (pop_of_interest, ''.join(pops_reference)))
            path2trained = opath.join(slim_path, swifr_trained_path)
            outfile = opath.join(slim_path, out_path, 'neutral', parameter_model_name+'_classified')
            if not opath.exists(outfile):
                formatted_txt = txt.format(**{
                    'path2trained' : path2trained,
                    'allstats_file' : allstats_file,
                    'pi_vals' : ' '.join([str(x) for x in pi_vec]),
                    'outfile' : outfile,
                    'log_file' : opath.join(bash_path,parameter_model_name+'_classify_'+pop_of_interest),
                    })
                fp = open(opath.join(bash_path, parameter_model_name+'_classify_'+pop_of_interest+'.sh'),'w')
                fp.write(formatted_txt)
                fp.close()
                os.system('sbatch '+opath.join(bash_path, parameter_model_name+'_classify_'+pop_of_interest+'.sh'))
            
            #os.system('swifr_test --path2trained '+path2trained+' --file '+allstats_file+
            #   ' --pi '+' '.join([str(x) for x in pi_vec])+' --outfile '+outfile)











