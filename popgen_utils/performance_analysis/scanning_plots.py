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
import math

try:
    import importlib.resources as ilresources
except ImportError:
    try:
        import importlib_resources as ilresources
    except ImportError:
        raise ImportError('Must install backport of importlib_resources if not using Python >= 3.7')


def make_plots(project_name, swifr_out_path, swifr_train_path, model_name_neutral, model_name_sweep, sweep_pos, pop_of_interest, sim_length, numbins=20, data_path=None):
    if data_path is None:
        data_path = config.params()['paths']['data']
    base_path = opath.join(data_path, project_name)
    slim_path = opath.join(base_path, 'slim')
    file = open(opath.join(slim_path, swifr_train_path, 'component_stats.txt'))
    stats = file.read()
    file.close()
    stats = stats.strip().splitlines()

    df = read_files(project_name, swifr_out_path, swifr_train_path, model_name_neutral, model_name_sweep, sweep_pos, pop_of_interest, sim_length, data_path)
    for stat in stats:
        figure_outpath = os.path.join(slim_path, swifr_out_path)
        figure_title = stat+'_peakplot'
        make_peakplot(df, stat, figure_outpath, figure_title, sim_length, sweep_pos, numbins)



def read_files(project_name, swifr_out_path, swifr_train_path, model_name_neutral, model_name_sweep, sweep_pos, pop_of_interest, sim_length, data_path=None):

    if data_path is None:
        data_path = config.params()['paths']['data']
    base_path = opath.join(data_path, project_name)
    slim_path = opath.join(base_path, 'slim')

    #file = open(opath.join(slim_path, swifr_train_path, 'component_stats.txt'))
    #stats = file.read()
    #file.close()
    #stats = stats.strip().splitlines()

    #neutral and sweep params
    yaml_file_neutral = open(opath.join(slim_path,f'{model_name_neutral}.yaml'))
    params_neutral = yaml.load(yaml_file_neutral)
    yaml_file_sweep = open(opath.join(slim_path,f'{model_name_sweep}.yaml'))
    params_sweep = yaml.load(yaml_file_sweep)

    # #read in all neutral files
    # df_list_neutral = []
    # for sim in range(int(params_neutral['sims'])):
    #     parameter_model_name = (f'{model_name_neutral}_sim-{sim}')
    #     classified_file = opath.join(slim_path, swifr_out_path, 'neutral', parameter_model_name+'_'+pop_of_interest+'_classified')
    #     df_list_neutral.append(pd.read_csv(classified_file, header=0, delim_whitespace=True, na_values='-998'))
    # neutral_df = pd.concat(df_list_neutral)

    #read in all sweep files
    #df_list_sweeppos = []
    #df_list_linked = []
    df_list = []
    for scoeff in params_sweep['selection_coefficient']:
        for time in params_sweep['sweep_time']:
            parameter_model_name = (f'{model_name_sweep}_coeff-{scoeff}_'
                                        f'pop-{pop_of_interest}_start-{time}')
            classified_file = opath.join(slim_path, swifr_out_path, 'sweep', parameter_model_name+'_classified')
            df = pd.read_csv(classified_file, header=0, delim_whitespace=True, na_values='-998')
            #sw = df.loc[df['pos'] == sweep_pos]
            #li = df.loc[df['pos'] != sweep_pos]
            #df_list_sweeppos.append(sw)
            #df_list_linked.append(li)
            df_list.append(df)
    #sweep_df = pd.concat(df_list_sweeppos)
    #linked_df = pd.concat(df_list_linked)
    df = pd.concat(df_list)

    return df

def make_peakplot(df, stat, figure_outpath, figure_title, sim_length, sweep_pos, numbins):
    #note: assumes sweep_pos is in the middle of the simulated region
    sweepvals = df.loc[df['pos'] == sweep_pos]
    sweepvals = sorted(sweepvals[stat].tolist())
    
    flankingdf = df.loc[df['pos'] != sweep_pos]
    #create a list of dataframes for flanking in 20 bins
    binlength = float(sim_length)/numbins
    flankingvals = []
    
    for i in range(numbins):
        lowerlimit = binlength*i
        upperlimit = binlength*(i+1)
        df_bin = df.loc[(df['pos']>lowerlimit) & (df['pos']<upperlimit)]
        flankingvals.append(sorted(df_bin[stat].tolist()))

    #look at 1st, 25th, 50th, 75th, 99th percentile
    percentile_1 = []
    percentile_25 = []
    percentile_50 = []
    percentile_75 = []
    percentile_99 = []


    #left flanking first:
    for i in range(numbins/2):
        percentile_1.append(flankingvals[i][math.floor(.01*length(flankingvals[i]))])
        percentile_25.append(flankingvals[i][math.floor(.25*length(flankingvals[i]))])
        percentile_50.append(flankingvals[i][math.floor(.5*length(flankingvals[i]))])
        percentile_75.append(flankingvals[i][math.floor(.75*length(flankingvals[i]))])
        percentile_99.append(flankingvals[i][math.floor(.99*length(flankingvals[i]))])
    #middle
    percentile_1.append(sweepvals[math.floor(.01*length(sweepvals))])
    percentile_25.append(sweepvals[math.floor(.25*length(sweepvals))])
    percentile_50.append(sweepvals[math.floor(.5*length(sweepvals))])
    percentile_75.append(sweepvals[math.floor(.75*length(sweepvals))])
    percentile_99.append(sweepvals[math.floor(.99*length(sweepvals))])
    #right flanking:
    for i in range(numbins/2, numbins):
        percentile_1.append(flankingvals[i][math.floor(.01*length(flankingvals[i]))])
        percentile_25.append(flankingvals[i][math.floor(.25*length(flankingvals[i]))])
        percentile_50.append(flankingvals[i][math.floor(.5*length(flankingvals[i]))])
        percentile_75.append(flankingvals[i][math.floor(.75*length(flankingvals[i]))])
        percentile_99.append(flankingvals[i][math.floor(.99*length(flankingvals[i]))])

    #values for the x-axis

    #left flanking
    for i in range(numbins/2):
        lowerlimit = binlength*i
        upperlimit = binlength*(i+1)
        x_axis_values.append((lowerlimit+upperlimit)/2)
    x_axis_values.append(sweep_pos)
    for i in range(numbins/2, numbins):
        lowerlimit = binlength*i
        upperlimit = binlength*(i+1)
        x_axis_values.appennd((lowerlimit+upperlimit)/2)

    plt.plot(x_axis_values, percentile_1, 'g-')
    plt.plot(x_axis_values, percentile_25, 'p-')
    plt.plot(x_axis_values, percentile_50, 'k-')
    plt.plot(x_axis_values, percentile_75, 'p-')
    plt.plot(x_axis_values, percentile_99, 'g-')
    plt.xlabel('Position in Simulated Genome Region')
    plt.ylabel(stat)
    plt.savefig(os.path.join(figure_outpath, figure_title+'.pdf'), bbox_inches='tight')
    plt.clf()




    



















