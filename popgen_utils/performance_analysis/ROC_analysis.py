import multiprocessing as mp
import subprocess
from matplotlib import pyplot as plt

from popgen_utils import config
from popgen_utils.misc import hashing

from datetime import datetime
import os
import os.path as opath
import yaml
import pandas as pd
try:
    import importlib.resources as ilresources
except ImportError:
    try:
        import importlib_resources as ilresources
    except ImportError:
        raise ImportError('Must install backport of importlib_resources if not using Python >= 3.7')

#read in classified examples, return dataframes, one for sweep simulations, one for neutral
# def read_classified_file(filename, classes):
#   df = pd.read_csv(filename, header=0, delim_whitespace=True)
#   names2use = ['pos']
#   for cl in classes:
#       names2use.append('P('+cl+')')
#   df = df[names2use]
#   return df

def read_classified_files_all(project_name, swifr_out_path, swifr_train_path, model_name_neutral, model_name_sweep, sweep_pos, pop_of_interest, data_path=None):
    '''
    project_name (str): name of project #gives the path
    swifr_out_path (str): where classified files live (top directories should correspond to neutral and sweep sims; relative to shared data path)
    swifr_train_path (str): where swifr was trained (relative to shared data path)
    model_name_neutral (str): name of yaml file for neutral files
    model_name_sweep (str): name of yaml file for sweep files
    sweep_pos (int) : position (1-indexed) of adaptive mutation
    pop_of_interest (str) : p1 e.g.
    '''

    #pops = ['p1','p2','p3']
    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path_neutral = opath.join(slim_path, model_name_neutral)
        slim_model_path_sweep = opath.join(slim_path,model_name_sweep)
        #bash_path = opath.join(slim_model_path,'bash')

        #extract classes in order -- NECESSARY?
        file = open(opath.join(slim_path,swifr_train_path,'classes.txt'))
        classes = file.read()
        file.close()
        classes = classes.strip().splitlines()

        #extract component stats -- NECESSARY?
        file = open(opath.join(slim_path, swifr_train_path, 'component_stats.txt'))
        stats = file.read()
        file.close()
        stats = stats.strip().splitlines()


        #neutral and sweep params
        yaml_file_neutral = open(opath.join(slim_path,f'{model_name_neutral}.yaml'))
        params_neutral = yaml.load(yaml_file_neutral)
        yaml_file_sweep = open(opath.join(slim_path,f'{model_name_sweep}.yaml'))
        params_sweep = yaml.load(yaml_file_sweep)

        #read in all neutral files
        df_list_neutral = []
        for sim in range(int(params_neutral['sims'])):
            parameter_model_name = (f'{model_name_neutral}_sim-{sim}')
            classified_file = opath.join(slim_path, swifr_out_path, 'neutral', parameter_model_name+'_'+pop_of_interest+'_classified')
            df_list_neutral.append(pd.read_csv(classified_file, header=0, delim_whitespace=True))
        neutral_df = pd.concat(df_list_neutral)

        #read in all sweep files
        df_list_sweeppos = []
        df_list_linked = []
        for scoeff in params_sweep['selection_coefficient']:
            for time in params_sweep['sweep_time']:
                parameter_model_name = (f'{model_name_sweep}_coeff-{scoeff}_'
                                            f'pop-{pop_of_interest}_start-{time}')
                classified_file = opath.join(slim_path, swifr_out_path, 'sweep', parameter_model_name+'_classified')
                df = pd.read_csv(classified_file, header=0, delim_whitespace=True)
                sw = df.loc[df['pos'] == sweep_pos]
                li = df.loc[df['pos'] != sweep_pos]
                df_list_sweeppos.append(sw)
                df_list_linked.append(li)
        sweep_df = pd.concat(df_list_sweeppos)
        linked_df = pd.concat(df_list_linked)

        return neutral_df, sweep_df, linked_df


def get_score_thresholds(list_of_scores):
    '''
    list_of_scores (list of floats, unsorted)

    returns: vector of 100 score thresholds (percentiles)
    '''
    scores = sorted(list_of_scores)
    indices = [int((int(x)*len(scores)/100)) for x in range(100)]
    #print(indices)
    return [scores[index] for index in indices]



def make_ROC_curves(project_name, swifr_out_path, swifr_train_path, model_name_neutral, model_name_sweep, sweep_pos, pop_of_interest, data_path=None):
    ''' 
    project_name (str): name of project #gives the path
    swifr_out_path (str): where classified files live (top directories should correspond to neutral and sweep sims; relative to shared data path)
    swifr_train_path (str): where swifr was trained (relative to shared data path)
    model_name_neutral (str): name of yaml file for neutral files
    model_name_sweep (str): name of yaml file for sweep files
    sweep_pos (int) : position (1-indexed) of adaptive mutation
    pop_of_interest (str) : p1 e.g.
    '''

    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')

    out_path = opath.join(slim_path,swifr_out_path)

    #read_classified_files_all
    [neutral_df, sweep_df, linked_df] = read_classified_files_all(project_name, 
        swifr_out_path, swifr_train_path, model_name_neutral, model_name_sweep, sweep_pos, pop_of_interest)
        
    file = open(opath.join(slim_path, swifr_train_path, 'component_stats.txt'))
    stats = file.read()
    file.close()
    stats = stats.strip().splitlines() 



    #for each statistic (read in from swifr_train_path), collect scores for neutral/linked/sweep --> get_score_thresholds
    stat2thresholds = {}
    for stat in stats:
        neutral_scores = neutral_df[stat].tolist()
        sweep_scores = sweep_df[stat].tolist()
        linked_scores = linked_df[stat].tolist()

        stat2thresholds[stat] = get_score_thresholds(neutral_scores+sweep_scores+linked_scores)
    neutral_psweep_aode = neutral_df['P(sweep)'].tolist()
    sweep_psweep_aode = sweep_df['P(sweep)'].tolist()
    linked_psweep_aode = linked_df['P(sweep)'].tolist()
    neutral_pliinked_aode = neutral_df['P(linked)'].tolist()
    linked_plinked_aode = neutral_df['P(linked)'].tolist()
    stat2thresholds['P(sweep)'] = get_score_thresholds(neutral_psweep_aode+sweep_psweep_aode+linked_psweep_aode)
    stat2thresholds['P(linked)'] = get_score_thresholds(neutral_pliinked_aode+linked_plinked_aode)


    #make ROC for sweep vs neutral
    stat2rates = {stat:[[],[]] for stat in stats}
    for stat in stats:
        [tp_rates, fp_rates] = get_tprate_fprate(neutral_df, sweep_df, stat, stat2thresholds[stat])
        stat2rates[stat][0] = tp_rates
        stat2rates[stat][1] = fp_rates
    [aode_tprates, aode_fprates] = get_tprate_fprate(neutral_df, sweep_df, 'P(sweep)', stat2thresholds['P(sweep)'])
    stat2rates['AODE'] = [aode_tprates, aode_fprates]

    plot_ROC(stat2rates, out_path, title='Sweep v Neutral')

    #make ROC for sweep vs linked
    stat2rates = {stat:[[],[]] for stat in stats}
    for stat in stats:
        [tp_rates, fp_rates] = get_tprate_fprate(linked_df, sweep_df, stat, stat2thresholds[stat])
        stat2rates[stat][0] = tp_rates
        stat2rates[stat][1] = fp_rates
    [aode_tprates, aode_fprates] = get_tprate_fprate(linked_df, sweep_df, 'P(sweep)', stat2thresholds['P(sweep)'])
    stat2rates['AODE'] = [aode_tprates, aode_fprates]

    plot_ROC(stat2rates, out_path, title='Sweep v Linked')

    #make ROC for linked vs neutral
    stat2rates = {stat:[[],[]] for stat in stats}
    for stat in stats:
        [tp_rates, fp_rates] = get_tprate_fprate(neutral_df, linked_df, stat, stat2thresholds[stat])
        stat2rates[stat][0] = tp_rates
        stat2rates[stat][1] = fp_rates
    [aode_tprates, aode_fprates] = get_tprate_fprate(neutral_df, linked_df, 'P(linked)', stat2thresholds['P(linked)'])
    stat2rates['AODE'] = [aode_tprates, aode_fprates]

    plot_ROC(stat2rates, out_path, title='Linked v Neutral')

def get_tprate_fprate(dataframe_neg, dataframe_pos, stat, thresholds, negate=False):
    '''
    dataframe_neg: dataframe for "negatives" 
    dataframe_pos: dataframe for "positives"
    stat: statisitc of interest
    thresholds: statistic-specific thresholds
    '''
    fps = [0 for i in range(len(thresholds))]
    tps = [0 for i in range(len(thresholds))]
    fns = [0 for i in range(len(thresholds))]
    tns = [0 for i in range(len(thresholds))]
    tp_rates = [0 for i in range(len(thresholds))]
    fp_rates = [0 for i in range(len(thresholds))]

    for i in range(len(thresholds)):
        thresh = thresholds[i]
        fps[i] = len(dataframe_neg[dataframe_neg[stat]>thresh])
        tns[i] = len(dataframe_neg[dataframe_neg[stat]<=thresh])
        fns[i] = len(dataframe_pos[dataframe_pos[stat]<=thresh])
        tps[i] = len(dataframe_pos[dataframe_pos[stat]>thresh])

        tp_rates[i] = float(tps[i])/(tps[i]+fns[i])
        fp_rates[i] = float(fps[i])/(fps[i]+tns[i])

    return [tp_rates, fp_rates]

def plot_ROC(stat2rates, out_path, title):
    '''
    stat2rates: dict mapping statistics (including 'AODE') to [tp_rates, fp_rates] #note: as of now, can only have 8 stats
    title (str): title for plot (and name for saving)
    out_path (str): directory to save plot in (then use title for file name)
    '''
    colors2use = ['#1b9e77','#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']
    color_index = -1
    statlist = []
    for stat in stat2rates.keys():
        statlist.append(stat)
        color_index += 1
        plt.plot(stat2rates[stat][1], stat2rates[stat][0], 'o-', color=colors2use[color_index], linewidth=2)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
    plt.legend(statlist)
    plt.savefig(opath.join(out_path, title+'.pdf'), format='pdf')
    plt.clf()




