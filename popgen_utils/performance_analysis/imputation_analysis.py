import multiprocessing as mp
import subprocess
from matplotlib import pyplot as plt

from popgen_utils import config
from popgen_utils.misc import hashing

from datetime import datetime
import numpy as np
import os
import os.path as opath
import yaml
import pandas as pd
import scipy.stats as ss
import random
try:
    import importlib.resources as ilresources
except ImportError:
    try:
        import importlib_resources as ilresources
    except ImportError:
        raise ImportError('Must install backport of importlib_resources if not using Python >= 3.7')


#read in neutral files: for each triplet of consecutive positions with ihs defined (or xpehh), extract middle value, and mean of left and right. 
#do the same for sweep files
#plot the two, calculate the pearson correlation

def read_files_neutral(project_name, model_name, pop_of_interest, figure_out_path, data_path=None):
    '''
    figure_out_path (str): where figure goes, relative to slim_path
    '''
    pops = ['p1','p2','p3']
    refpops = [pop for pop in pops if pop != pop_of_interest]
    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path,model_name)

        yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
        params = yaml.load(yaml_file)

        imputed_values = {}
        actual_values = {}

        for sim in range(int(params['sims'])):
            parameter_model_name = (f'{model_name}_sim-{sim}')
            allstats_file = opath.join(slim_model_path, parameter_model_name+'_'+pop_of_interest+'_'+''.join(refpops)+'_allstats.txt')
            df = pd.read_csv(allstats_file, header=0, delim_whitespace=True, na_values='-998')
            stats = df.columns[2:]
            for stat in stats:
                vals = df[stat].tolist()
                for i in range(len(vals)-3):
                    if np.isnan(vals[i])==False and np.isnan(vals[i+1])==False and np.isnan(vals[i+2])==False:
                        if random.random()<0.01:
                            actual_value = vals[i+1]
                            imputed_value = float(vals[i]+vals[i+2])/2
                            if stat in imputed_values.keys():
                                actual_values[stat].append(actual_value)
                                imputed_values[stat].append(imputed_value)
                            else:
                                actual_values[stat] = [actual_value]
                                imputed_values[stat] = [imputed_value]
       
        for stat in imputed_values.keys():

            plt.plot(actual_values[stat], imputed_values[stat], 'o')
            plt.xlabel('Actual Statistic Value')
            plt.ylabel('Imputed Statistic Value')
            plt.title(stat)
            plt.savefig(opath.join(slim_path, figure_out_path,'imputation_'+stat+'_neutral.pdf'))
            plt.clf()
            print(stat+' correlation: '+str(ss.pearsonr(actual_values[stat], imputed_values[stat])[0]))


def read_files_sweep(project_name, model_name, pop_of_interest, figure_out_path, data_path=None):

    pops = ['p1','p2','p3']
    refpops = [pop for pop in pops if pop != pop_of_interest]

    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path,model_name)

        yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
        params = yaml.load(yaml_file)

        imputed_values = {}
        actual_values = {}

        for scoeff in params['selection_coefficient']:
            for time in params['sweep_time']:
                parameter_model_name = (f'{model_name}_coeff-{scoeff}_'
                                            f'pop-{pop_of_interest}_start-{time}')
                allstats_file = opath.join(slim_model_path, parameter_model_name+'_'+pop_of_interest+'_'+''.join(refpops)+'_allstats.txt')

                df = pd.read_csv(allstats_file, header=0, delim_whitespace=True, na_values='-998')
                stats = df.columns[2:]
                for stat in stats:
                    vals = df[stat].tolist()
                    for i in range(len(vals)-3):
                        if np.isnan(vals[i])==False and np.isnan(vals[i+1])==False and np.isnan(vals[i+2])==False:
                            if random.random()<0.01:
                                actual_value = vals[i+1]
                                imputed_value = float(vals[i]+vals[i+2])/2
                                if stat == 'ihs':
                                    actual_value = abs(vals[i+1])
                                    imputed_value = float(abs(vals[i]+vals[i+2])/2)
                                if stat in imputed_values.keys():
                                    actual_values[stat].append(actual_value)
                                    imputed_values[stat].append(imputed_value)
                                else:
                                    actual_values[stat] = [actual_value]
                                    imputed_values[stat] = [imputed_value]
           
            for stat in imputed_values.keys():

                plt.plot(actual_values[stat], imputed_values[stat], 'o')
                plt.xlabel('Actual Statistic Value')
                plt.ylabel('Imputed Statistic Value')
                plt.title(stat)
                plt.savefig(opath.join(slim_path, figure_out_path,'imputation_'+stat+'_sweep.pdf'))
                plt.clf()
                print(stat+' correlation: '+str(ss.pearsonr(actual_values[stat], imputed_values[stat])))

