import multiprocessing as mp
import subprocess
import pandas as pd

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

#need to use .hap and .map files: .map files give the locations (for windowing), hap has the 0s and 1s

def run_window_stats(project_name, model_name, type, windowsize, data_path=None):
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

    if type == 'sweep':
        scoeffs = params['selection_coefficient']
        times = params['sweep_time']
        pops_of_interest = params['sweep_population']
        for coeff in scoeffs:
            for time in times:
                for pop in pops_of_interest:
                    parameter_model_name = (f'{model_name}_coeff-{coeff}_'f'pop-{pop}_start-{time}')
                    filename = opath.join(slim_model_path, f'{parameter_model_name}')
                    outpath = filename+'_window_allstats.txt'
                    df = compute_sfs(filename, pop, windowsize)
                    df.to_csv(outpath, sep=' ')

    elif type == 'neutral':
        sims = params['sims']
        for sim in range(int(sims)):
            for pop in ['1','2','3']:
                parameter_model_name = (f'{model_name}_sim-{sim}')
                filename = opath.join(slim_model_path, parameter_model_name)
                outpath = filename+'_window_allstats.txt'
                df = compute_sfs(filename, pop, windowsize)
                df.to_csv(outpath, sep=' ')





def compute_sfs(filename, pop, windowsize):

    windowsize=int(windowsize)
    mapfile = filename+'_map.txt'
    hapfile = filename+'_p'+pop+'.hap' #only hapfile for pop of interest
    df1 = pd.read_csv(mapfile, sep=' ', header=None, names=['chrom','snpname','mapdist','pos'])
    positions = df1['pos'].tolist()
    df2 = pd.read_csv(hapfile, sep=' ', header=None)
    nchroms = len(df2)
    dafs = df2.sum(axis=0).tolist() # in order, derived allele frequencies
    df1['dafs'] = dafs
    #go by windowsize:
    numwindows = float(positions[-1])/windowsize
    numwindows = int(numwindows) #floor

    Start = []
    End = []
    TW = []
    TPi = []
    TH = []
    TL = []
    TajD = []
    FWH = []
    ZE = []

    for window_i in range(numwindows):
        start_pos = window_i*windowsize
        end_pos = (window_i+1)*windowsize
        sub_df = df1[df1['pos'].isin(range(start_pos, end_pos))]
        dafs = sub_df['dafs'].tolist()
        etas = [len([x for x in dafs if int(x)==i]) for i in range(nchroms+1)] #note: includes fixed mutations

        Start.append(start_pos)
        End.append(end_pos)
        TW.append(thetaW(etas, nchroms))
        TPi.append(thetaPi(etas, nchroms))
        TH.append(thetaH(etas, nchroms))
        TL.append(thetaL(etas, nchroms))
        TajD.append(TajimaD(etas, nchroms))
        FWH.append(FayWuH(etas, nchroms))
        ZE.append(ZengE(etas, nchroms))

    df = pd.DataFrame(data={'Start':Start, 'End':End, 'thetaW':TW, 'thetaPi':TPi, 'thetaH':TH, 'thetaL':TL, 
                             'TajimaD':TajD, 'FayWuH':FWH, 'ZengE':ZE})
    return df


def thetaW(etas, nchroms):
    '''
    :return: theta_W
    '''
    S = sum(etas[1:nchroms]) #only segregating sites
    a = sum([float(1)/i for i in range(1, nchroms)])
    return float(S)/a

def thetaPi(etas, nchroms): 
    '''
    :return: theta_pi
    '''

    a = float(2)/(nchroms*(nchroms-1))
    S = 0
    for i in range(1,nchroms):
        S += i*(nchroms-i)*etas[i]
    return a*S

def thetaH(etas, nchroms): 
    '''
    :return: theta_H (Fay & Wu)
    '''

    a = float(2)/(nchroms*(nchroms-1))
    S = 0
    for i in range(1,nchroms):
        S += i*i*etas[i]
    return a*S

def thetaL(etas, nchroms):
    '''
    :return: theta_L (Zeng)
    '''
    a = float(1)/(nchroms-1)
    S = 0
    for i in range(1,nchroms):
        S += i*etas[i]
    return a*S

def TajimaD(etas, nchroms):
    '''
    :return: Tajima's D, unnormalized (numerator only)
    '''
    return thetaPi(etas, nchroms)-thetaW(etas, nchroms)

def FayWuH(etas, nchroms):
    '''
    :return: Fay & Wu's H, unnormalized (numerator only)
    '''
    return thetaPi(etas, nchroms)-thetaH(etas, nchroms)

def ZengE(etas, nchroms):
    '''
    :return: Zeng's E, unnormalized (numerator only)
    '''
    return thetaL(etas, nchroms)-thetaW(etas, nchroms)  



