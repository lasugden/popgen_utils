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
#from popgen_utils.hmm import oscar_scripts
from popgen_utils.hmm import ode_hmm
from popgen_utils.hmm.ode_hmm import Stats, Mix1D


def run_ode_hmm_on_training_data(project_name, model_name, type, pop_of_interest,
    swifr_trained_path, out_path, data_path=None):
    '''
    project_name (str): name of project (directory)
    model_name (str): e.g. gravel_neutral_400kb
    type (str): 'sweep' or 'neutral'
    pop_of_interest (str): in the format 'pN' where N is integer in 1-3
    swifr_trained_path (str): path where trained model lives
    out_path (str): path for the output (relative path from common directory)
    '''

    #make out path if it doesn't exist

    pops = ['p1','p2','p3']
    if data_path is None:
        data_path = config.params()['paths']['data']
    base_path = opath.join(data_path, project_name)
    slim_path = opath.join(base_path, 'slim')
    slim_model_path = opath.join(slim_path, model_name)
    bash_path = opath.join(slim_model_path,'bash')
    output_path = opath.join(slim_path, out_path)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)

    path2trained = opath.join(slim_path, swifr_trained_path)
    HMM = ode_hmm.ODE_HMM(path2trained) #class ODE_HMM -- can use functions within
    print('mkdir -p '+output_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)


    if type == 'sweep':
        for scoeff in params['selection_coefficient']:
            for time in params['sweep_time']:
                pops_reference = [x for x in pops if x != pop_of_interest]
                parameter_model_name = (f'{model_name}_coeff-{scoeff}_'
                                            f'pop-{pop_of_interest}_start-{time}')
                allstats_file = opath.join(slim_model_path, 
                    parameter_model_name+'_%s_%s_allstats.txt' % (pop_of_interest, ''.join(pops_reference)))
                [positions, Svec] = HMM.read_Svec(allstats_file)
                HMM.viterbi(Svec, positions, output_path, parameter_model_name, plotpaths=True, classify=True) 
                HMM.many_backtraces(Svec, 1000, positions, output_path, parameter_model_name, plotpaths=True, plotdensity=True, classify=True, plot_classify=True)

    elif type == 'neutral':
        for sim in range(int(params['sims'])):
            pops_reference = [x for x in pops if x!= pop_of_interest]
            parameter_model_name = (f'{model_name}_sim-{sim}')
            allstats_file = opath.join(slim_model_path, 
                parameter_model_name+'_%s_%s_allstats.txt' % (pop_of_interest, ''.join(pops_reference)))
            [positions, Svec] = HMM.read_Svec(allstats_file)
            HMM.viterbi(Svec, positions, output_path, parameter_model_name, plotpaths=True, classify=True)
            HMM.many_backtraces(Svec, 1000, positions, output_path, parameter_model_name, plotpaths=True, plotdensity=True, classify=True, plot_classify=True)





