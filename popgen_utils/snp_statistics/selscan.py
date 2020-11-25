import multiprocessing as mp
import subprocess

from datetime import datetime
import os
import os.path as opath
import yaml
import random
try:
    import importlib.resources as ilresources
except ImportError:
    try:
        import importlib_resources as ilresources
    except ImportError:
        raise ImportError('Must install backport of importlib_resources if not using Python >= 3.7')

from popgen_utils import config
from popgen_utils.misc import hashing

from popgen_utils.snp_statistics import oscar_scripts

def normalize_ihs(project_name, model_name, data_path=None):
    # first normalize all neutral simulations together
    pops = ['p1','p2','p3']
    
    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path, model_name)
        bash_path = opath.join(slim_model_path,'bash')
        if not opath.exists(bash_path):
            os.mkdir(bash_path)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)    

    txt_neutral = ilresources.open_text(oscar_scripts, 'norm_neutral.txt').read()
    txt_sweep = ilresources.open_text(oscar_scripts, 'norm_sweep.txt').read()

    sims = params['sims']
    for pop in pops:
        neutral_files = ''
        for sim in range(int(sims)):
            parameter_model_name = (f'{model_name}_sim-{sim}')
            filename = opath.join(slim_model_path,f'{parameter_model_name}_{pop}.ihs.out')
            neutral_files = neutral_files+' '+filename
            formatted_txt = txt_neutral.format(**{
                'file_list' : neutral_files,
                'log_file' : 'ihs_norm_neutral_'+pop
                })
        fp = open(opath.join(bash_path,'normalize_ihs_neutral_'+pop+'.sh'),'w')
        fp.write(formatted_txt)
        fp.close()
        os.system('sbatch '+opath.join(bash_path,'normalize_ihs_neutral.sh'))

        #for this population of interest, normalize sweep files


        for coeff in params['selection_coefficient']:
            for time in params['sweep_time']:
            #for pop in params['sweep_population']:
                parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                            f'pop-{pop}_start-{time}')
                filename = opath.join(slim_model_path,f'{parameter_model_name}_{pop}.ihs.out')
                if not os.path.isfile(opath.join(slim_model_path,f'{parameter_model_name}_{pop}.ihs.out.100bins.norm')):
                    formatted_txt = txt_sweep.format(**{
                        'sweep_file' : filename,
                        'neutral_file_list' : neutral_files,
                        'log_file' : 'ihs_norm_sweep_'+parameter_model_name
                        })
                    fp = open(opath.join(bash_path,'normalize_ihs_sweep_'+parameter_model_name+'.sh'),'w')
                    fp.write(formatted_txt)
                    fp.close()
                    os.system('sbatch '+opath.join(bash_path,'normalize_ihs_sweep_'+parameter_model_name+'.sh'))



def run_ihs(project_name, model_name, type, data_path=None):

    pops = ['p1','p2','p3']

    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path, model_name)
        bash_path = opath.join(slim_model_path,'bash')
        if not opath.exists(bash_path):
            os.mkdir(bash_path)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)    

    txt = ilresources.open_text(oscar_scripts, 'run_ihs.txt').read()

    if type == 'sweep':
        scoeffs = params['selection_coefficient']
        times = params['sweep_time']
        pops_of_interest = params['sweep_population']
        for coeff in scoeffs:
            for time in times:
                for pop in pops_of_interest:
                    parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                            f'pop-{pop}_start-{time}')
                    if not opath.isfile(opath.join(slim_model_path,f'{parameter_model_name}_{pop}.ihs.out')) or os.path.getsize(opath.join(slim_model_path,f'{parameter_model_name}_{pop}.ihs.out'))<2:

                        formatted_txt = txt.format(**{
                            'hap_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}.hap'),
                            'map_file' : opath.join(slim_model_path,f'{parameter_model_name}_map.txt'),
                            'out_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}'),
                            'log_file' : opath.join(bash_path,f'{parameter_model_name}_{pop}'),
                            })
                        fp = open(opath.join(bash_path,f'{parameter_model_name}_{pop}.sh'),'w')
                        fp.write(formatted_txt)
                        fp.close()
                        os.system('sbatch '+opath.join(bash_path,f'{parameter_model_name}_{pop}.sh'))

    elif type == 'neutral':
        sims = params['sims']
        for sim in range(int(sims)):
            for pop in pops:
                parameter_model_name = (f'{model_name}_sim-{sim}')
                if not opath.isfile(opath.join(slim_model_path,f'{parameter_model_name}_{pop}.ihs.out')):
                    formatted_txt = txt.format(**{
                        'hap_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}.hap'),
                        'map_file' : opath.join(slim_model_path,f'{parameter_model_name}_map.txt'),
                        'out_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}'),
                        'log_file' : opath.join(bash_path,f'{parameter_model_name}_{pop}'),
                        })
                    fp = open(opath.join(bash_path,f'{parameter_model_name}_{pop}.sh'),'w')
                    fp.write(formatted_txt)
                    fp.close()
                    os.system('sbatch '+opath.join(bash_path,f'{parameter_model_name}_{pop}.sh'))



def run_xpehh(project_name, model_name, type, data_path=None):

    pops = ['p1','p2','p3']

    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path, model_name)
        bash_path = opath.join(slim_model_path,'bash')
        if not opath.exists(bash_path):
            os.mkdir(bash_path)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)    

    txt = ilresources.open_text(oscar_scripts, 'run_xpehh.txt').read()

    if type == 'sweep':
        scoeffs = params['selection_coefficient']
        times = params['sweep_time']
        pops_of_interest = params['sweep_population']
        for coeff in scoeffs:
            for time in times:
                for pop in pops_of_interest:
                    for refpop in [x for x in pops if x!= pop]:
                        parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                                f'pop-{pop}_start-{time}')
                        if not opath.isfile(opath.join(slim_model_path,f'{parameter_model_name}_{pop}_{refpop}.xpehh.out')):

                            formatted_txt = txt.format(**{
                                'hap_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}.hap'),
                                'ref_file' : opath.join(slim_model_path,f'{parameter_model_name}_{refpop}.hap'),
                                'map_file' : opath.join(slim_model_path,f'{parameter_model_name}_map.txt'),
                                'out_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}_{refpop}'),
                                'log_file' : opath.join(bash_path,f'{parameter_model_name}_{pop}_{refpop}'),
                                })
                            fp = open(opath.join(bash_path,f'{parameter_model_name}_{pop}_{refpop}.sh'),'w')
                            fp.write(formatted_txt)
                            fp.close()
                            os.system('sbatch '+opath.join(bash_path,f'{parameter_model_name}_{pop}_{refpop}.sh'))

    elif type == 'neutral':
        sims = params['sims']
        for sim in range(int(sims)):
            for pop in pops:
                for refpop in pops:
                    if pops.index(refpop) > pops.index(pop): #only run p1p2 e.g., not p2p1
                        parameter_model_name = (f'{model_name}_sim-{sim}')
                        if not opath.isfile(opath.join(slim_model_path,f'{parameter_model_name}_{pop}_{refpop}.xpehh.out')):

                            formatted_txt = txt.format(**{
                                'hap_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}.hap'),
                                'ref_file' : opath.join(slim_model_path,f'{parameter_model_name}_{refpop}.hap'),
                                'map_file' : opath.join(slim_model_path,f'{parameter_model_name}_map.txt'),
                                'out_file' : opath.join(slim_model_path,f'{parameter_model_name}_{pop}_{refpop}'),
                                'log_file' : opath.join(bash_path,f'{parameter_model_name}_{pop}_{refpop}'),
                                })
                            fp = open(opath.join(bash_path,f'{parameter_model_name}_{pop}_{refpop}.sh'),'w')
                            fp.write(formatted_txt)
                            fp.close()
                            os.system('sbatch '+opath.join(bash_path,f'{parameter_model_name}_{pop}_{refpop}.sh'))









