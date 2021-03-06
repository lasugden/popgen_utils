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

import argparse
import pandas as pd


def read_file_mapfile(directory_path, parameter_model_name, pop_of_interest=None, refpop=None):
    mappath = os.path.join(directory_path, parameter_model_name+'_map.txt')
    if not os.path.exists(mappath):
        return None
    df = pd.read_csv(mappath, skiprows=0, header=None, delim_whitespace=True,
            names=['chrom_name', 'locus_name', 'mappos', 'pos'])
    return df[['locus_name','pos']]

def read_file_ihs(directory_path, parameter_model_name, pop_of_interest, refpop=None, key_column='iHS'):
    """Read ihs files

    Args:
        directory_path (str): path to directory containing files
        seed (int): seed used to create the files
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        refpop (str, unused): in the format 'pN' where N is an integer in 1-3
        key_column (str, optional): Name of the key column used by the file. Defaults to 'iHS'.
    """
    ihs_path = os.path.join(directory_path, parameter_model_name+'_%s.ihs.out.100bins.norm' % (pop_of_interest))
    if not os.path.exists(ihs_path):
        return pd.DataFrame(columns=['pos', key_column])
    if os.stat(os.path.join(directory_path, parameter_model_name+'_%s.ihs.out.100bins.norm' % (pop_of_interest))).st_size == 0:
        return pd.DataFrame(columns=['pos', key_column])
    df = pd.read_csv(ihs_path, skiprows=0, header=None, delim_whitespace=True, usecols=range(7),
                     names=['locus_name', 'pos', 'freq', 'ihh1', 'ihh2', 'ihs_unnormalized', key_column])
    return df[['pos', key_column]]


def read_file_xpehh(directory_path, parameter_model_name, pop_of_interest, refpop, key_column='XP-EHH'):
    """Read XPEHH files

    Args:
        directory_path (str): path to directory containing files
        seed (int): seed used to create the files
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        refpop (str): in the format 'pN' where N is an integer in 1-3
        key_column (str, optional): Name of the key column used by the file. Defaults to 'iHS'.
    """
    xpehh_path = os.path.join(directory_path, parameter_model_name+'_%s_%s.xpehh.out' % (pop_of_interest, refpop))
    xpehh_path_alt = os.path.join(directory_path, parameter_model_name+'_%s_%s.xpehh.out' % (refpop, pop_of_interest))
    if not os.path.exists(xpehh_path) and not os.path.exists(xpehh_path_alt):
        return pd.DataFrame(columns=['pos', key_column])

    if os.path.exists(xpehh_path):
        df = pd.read_csv(xpehh_path, skiprows=1, header=None, delim_whitespace=True,
                     names=['locus_name', 'pos', 'gpos', 'p1', 'ihh1', 'p2', 'ihh2', key_column])
    else:
        df = pd.read_csv(xpehh_path_alt, skiprows=1, header=None, delim_whitespace=True,
                         names=['locus_name', 'pos', 'gpos', 'p1', 'ihh1', 'p2', 'ihh2', key_column])
        df[key_column] = -1*df[key_column]
    return df[['pos', key_column]]


def read_file_isafe(directory_path, parameter_model_name, pop_of_interest, refpop=None, key_column='iSAFE'):
    """Read isafe files

    Args:
        directory_path (str): path to directory containing files
        seed (int): seed used to create the files
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        pops_reference (str, unused): in the format 'pN' where N is an integer in 1-3
        key_column (str, optional): Name of the key column used by the file. Defaults to 'iHS'.
    """
    isafe_path = os.path.join(directory_path, parameter_model_name+'_%s.iSAFE.out' % (pop_of_interest))
    if not os.path.exists(isafe_path):
        return pd.DataFrame(columns=['pos', key_column])
    df = pd.read_csv(isafe_path, skiprows=1, header=None, delim_whitespace=True,
                     names=['pos', key_column, 'daf'])
    # Match locations for 0-indexing and 1-indexing
    df['pos'] -= 1
    return df[['pos', key_column]]


def read_file_fst(directory_path, parameter_model_name, pop_of_interest, refpop, key_column='Fst'):
    """Read fst files

    Args:
        directory_path (str): path to directory containing files
        seed (int): seed used to create the files
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        refpop (str): in the format 'pN' where N is an integer in 1-3
        key_column (str, optional): Name of the key column used by the file. Defaults to 'iHS'.
    """
    #fst_path = os.path.join(directory_path, '%i_%s%s.weir.fst' % (seed, pop_of_interest, refpop))
    fst_path = os.path.join(directory_path, parameter_model_name+'_%s_%s.weir.fst' % (pop_of_interest, refpop))
    fst_path_alt = os.path.join(directory_path, parameter_model_name+'_%s_%s.weir.fst' % (refpop, pop_of_interest))
    if not os.path.exists(fst_path) and not os.path.exists(fst_path_alt):
        return pd.DataFrame(columns=['pos', key_column])
    if os.path.exists(fst_path):
        df = pd.read_csv(fst_path, skiprows=1, header=None, delim_whitespace=True,
                     names=['chrom', 'pos', key_column])
    else:
        df = pd.read_csv(fst_path_alt, skiprows=1, header=None, delim_whitespace=True,
                     names=['chrom', 'pos', key_column])
    # Match locations for 0-indexing and 1-indexing
    df['pos'] -= 1
    return df[['pos', key_column]]


# def get_seeds(path):
#     """
#     From a directory, get all seed values by searching over names

#     Args:
#         path (str): path to search

#     Returns:
#         list of ints: the seeds that can be iterated over

#     """
#     seeds = os.listdir(path)
#     seeds = list(set([int(x[:x.find('_')]) for x in seeds if x.find('100bins.norm') > -1]))
#     return seeds


def read_files(directory_path, parameter_model_name, pop_of_interest, pops_reference):
    """
    Read in IHS, XPEHH, ISAFE, and FST files, merge on position, and return all statistics

    Args:
        directory_path (str): path to directory containing files
        seed (int): seed used to create the files
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        pops_reference (str/list of str): in the format 'pN' where N is an integer in 1-3

    Returns:
        Pandas DataFrame: the merged output of the statistics

    """
    df = read_file_mapfile(directory_path, parameter_model_name)
    if df is None:
        raise NotImplementedError('Could not find map file. We do not yet account for naming of loci without reading in mapfile data.')

    #refpop_num = 0

    for refpop in pops_reference:
        #refpop_num += 1
        df_fst = read_file_fst(directory_path, parameter_model_name, pop_of_interest, refpop, key_column='fst_'+refpop)
        #if df_fst is None and refpop_num == 1:
        #    raise NotImplementedError('Could not find fst file. We do not yet account for naming of loci without reading in Fst data.')
        #if refpop_num == 1:
        #    df = df_fst
        if df_fst is not None:
            df = df.merge(df_fst, how='outer', on='pos')
        else:
            print('WARNING: FST file not found')    

    df_ihs = read_file_ihs(directory_path, parameter_model_name, pop_of_interest, refpop=None, key_column='ihs')
    if df_ihs is not None:
        df = df.merge(df_ihs, how='outer', on='pos')
        #raise NotImplementedError('Could not find ihs file. We do not yet account for naming of loci without reading in IHS data.')

    # Read in files

    for refpop in pops_reference:
        df_xpehh = read_file_xpehh(directory_path, parameter_model_name, pop_of_interest, refpop, key_column='xpehh_'+refpop )
        if df_xpehh is not None:
            df = df.merge(df_xpehh, how='outer', on='pos')
        else:
            print('WARNING: XPEHH file not found: refpop='+refpop)

    df_isafe = read_file_isafe(directory_path, parameter_model_name, pop_of_interest, refpop=None, key_column='isafe')
    if df_isafe is not None:
        df = df.merge(df_isafe, how='outer', on='pos')
    else:
        print('WARNING: iSAFE file not found')



    # Fill in NaNs appropriately
    #nan_names = df['locus_name'].isnull()
    #if nan_names.sum() > 0:
    #    df.loc[nan_names, 'locus_name'] = 'SNP_POS_' + df.loc[nan_names, 'pos'].apply(str)

    df = df.loc[df.isnull().sum(axis=1) < 4, :]
    df = df.fillna(-998)

    return df


def write_output(directory_path, parameter_model_name, pop_of_interest, pops_reference, df):
    """
    Write the merged output to a file

    Args:
        directory_path (str): directory into which the file should be saved
        seed (int): the seed used to create the file type OR combination of metadata
        pop_of_interest (str): of the type pN where N is 1+
        pops_reference (str/list of str): of the type pN where N is 1+
        df (Pandas DataFrame): the combined statistics for writing

    """
    #seed = '' if seed is None else seed
    file_name = parameter_model_name+'_%s_%s_allstats.txt' % (pop_of_interest, ''.join(pops_reference))
    file_path = os.path.join(directory_path, file_name)
    df.to_csv(file_path, sep='\t', index=False)


def write_allstats(input_directory,
                              parameter_model_name,
                              pop_of_interest,
                              pops_reference):
    """
    Iterate over all seeds in a directory and merge.

    Args:
        input_directory (str): the directory in which all of the input data lives
        output_directory (str): the directory to which the combined statistics should be written
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        pops_reference (str/list of str): in the format 'pN' where N is an integer in 1-3

    """
    pops_reference = [pops_reference] if isinstance(pops_reference, str) else pops_reference
    #seeds = get_seeds(input_directory)
    #for seed in seeds:

    df = read_files(input_directory, parameter_model_name, pop_of_interest, pops_reference)
    if df is None:
        print('WARNING: some files are missing for seed %i and population %s' % (seed, pop_of_interest))
    else:
        write_output(input_directory, parameter_model_name, pop_of_interest, pops_reference, df)


def write_allstats_single_snp(input_directory,
                                parameter_model_name,
                                pop_of_interest,
                                pops_reference,
                                sweep_pos,
                                additional_name_text=None):
    """
    Iterate over all seeds in a directory and merge.

    Args:
        input_directory (str): the directory in which all of the input data lives
        output_directory (str): the directory to which the combined statistics should be written
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        pops_reference (str/list of str): in the format 'pN' where N is an integer in 1-3
        sweep_pos (int): the position of the sweep to extract
        additional_name_text (str, optional): if set, this will be appended to filename

    """
    #pops_of_interest = [pops_of_interest] if isinstance(pops_of_interest, str) else pops_of_interest
    pops_reference = [pops_reference] if isinstance(pops_reference,str) else pops_reference
    #seeds = get_seeds(input_directory)
    combined = []
    #for population in pops_of_interest:
    #for refpop in pops_reference:
    #for seed in seeds:
    df = read_files(input_directory, parameter_model_name, pop_of_interest, pops_reference)
    if df is None:
        print('WARNING: some files are missing for seed %i and population %s' % (seed, pop_of_interest))
    else:
        combined.append(df.loc[df['pos'] == sweep_pos, :])

    if len(combined) > 0:
        df = pd.concat(combined, axis=0, ignore_index=True)
        write_output(input_directory, paramter_model_name+'_snp'+str(sweep_pos), pop_of_interest, pops_reference, df)


# def merge_seeds_over_stp(input_directory,
#                          output_directory,
#                          selection_strength,
#                          time_of_sweep,
#                          pop_of_sweep,
#                          pop_of_interest,
#                          pops_reference,
#                          sweep_pos=None):


def extract_compstats_region(project_name, model_name, pop_of_interest, sweep_pos_1, sweep_pos_2=None, additional_name_text=None, data_path=None):
    pops = ['p1','p2','p3']

    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path, model_name)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)

    scoeffs = params['selection_coefficient']
    times = params['sweep_time']
    #pops_of_interest = params['sweep_population']

    df_list = []
    #param_list = []

    for coeff in scoeffs:
        for time in times:
            #for pop in pops_of_interest:
            pops_reference = [x for x in pops if x!= pop_of_interest]
            parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                    f'pop-{pop_of_interest}_start-{time}')
            #param_list.append(parameter_model_name)
            allstats_file_name = parameter_model_name+'_%s_%s_allstats.txt' % (pop_of_interest, ''.join(pops_reference))
            df = pd.read_csv(opath.join(slim_model_path,allstats_file_name), header=0, delim_whitespace=True)
            if sweep_pos_2 == None:
                df_slice = df.loc[df['pos']==sweep_pos_1]
                df_slice['parameters'] = parameter_model_name
                df_list.append(df_slice)
            else:
                df_slice = df.loc[df['pos'].isin(range(sweep_pos_1, sweep_pos_2+1))]
                df_slice['parameters'] = parameter_model_name
                df_list.append(df_slice)
    region_allstats = pd.concat(df_list)
    #region_allstats['parameters'] = param_list

    file_name = pop_of_interest+'_'+additional_name_text+'_allstats.txt'
    file_path = os.path.join(slim_model_path, file_name)
    region_allstats.to_csv(file_path, sep='\t', index=False)




def merge_compstats(project_name, model_name, type, data_path=None, sweep_pos=None):
    """
    Iterate over populations, selection strengths, and sweep times and extract combined statistics.

    Args:
        input_directory (str): OUTER path that includes all simulations
        output_directory (str): path into which allstats will be written
        selection_strength (float/list of floats): selection strength
        time_of_sweep (int/list of ints): kiloyears of sweep
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        pop_of_sweep (str): the population in which the sweep occured
        pops_reference (str/list of str): in the format 'pN' where N is an integer in 1-3
        sweep_pos (int, optional): the position of the sweep to extract
            if None, do not extract a single sweep_pos, but rather save each file

    """
    pops = ['p1','p2','p3']
    if data_path is None:
        data_path = config.params()['paths']['data']
        base_path = opath.join(data_path, project_name)
        slim_path = opath.join(base_path, 'slim')
        slim_model_path = opath.join(slim_path, model_name)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)
 


    #selection_strength = [selection_strength] if not isinstance(selection_strength, list) else selection_strength
    #time_of_sweep = [time_of_sweep] if not isinstance(time_of_sweep, list) else time_of_sweep
    #pops_reference = [pops_reference] if isinstance(pops_reference,str) else pops_reference

    if type == 'sweep':
        scoeffs = params['selection_coefficient']
        times = params['sweep_time']
        pops_of_interest = params['sweep_population']
        for coeff in scoeffs:
            for time in times:
                for pop in pops_of_interest:
                    pops_reference = [x for x in pops if x!= pop]
                    parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                            f'pop-{pop}_start-{time}')
                    if sweep_pos is None:
                        write_allstats(slim_model_path, parameter_model_name, pop, pops_reference)
                    else:
                        write_allstats_single_snp(slim_model_path, parameter_model_name, pop, pops_reference, sweep_pos)

    elif type == 'neutral':
        sims = params['sims']
        for sim in range(int(sims)):
            for pop in pops:
                pops_reference = [x for x in pops if x!= pop]
                parameter_model_name = (f'{model_name}_sim-{sim}')
                write_allstats(slim_model_path, parameter_model_name, pop, pops_reference)

    # for sel in selection_strength:
    #     for t in time_of_sweep:
    #         sub_input = os.path.join(input_directory, pop_of_sweep, 'T%i' % t, 's%s' % sel)
    #         if sweep_pos is None:
    #             merge_all_seeds_and_write(sub_input, output_directory, pop_of_interest, pops_reference)
    #         else:
    #             merge_all_seeds_and_extract(sub_input,
    #                                         output_directory,
    #                                         pop_of_interest,
    #                                         pops_reference,
    #                                         sweep_pos,
    #                                         't%i_s%s' % (t, sel))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_path', help='Path to the directory to iterate over', type=str)
    parser.add_argument('--write_path', help='Path to the directory to write into. Must already exist.', type=str)
    parser.add_argument('--pop_of_interest', help='Comma-separated list of pops (pN; with no spaces).', type=str)
    parser.add_argument('--pops_reference', help='Comma-separated list of reference pops (pN; with no spaces).', type=str)
    parser.add_argument('--type', help='Can be neutral or sweep', type=str, choices=['neutral', 'sweep'])
    parser.add_argument('--times_of_sweep', help='Comma-separated list of integer times of sweep', type=str)
    parser.add_argument('--selection_strengths', help='Comma-separated list of float selection strengths', type=str)
    parser.add_argument('--pop_of_sweep', help='Population of sweep (pN; with no spaces)', type=str)
    parser.add_argument('--sweep_pos', help='Sweep position', default=500000, type=int)

    args = parser.parse_args()

    if args.type == 'neutral':
        merge_all_seeds_and_write(args.read_path,
                                  args.write_path,
                                  #args.pops_of_interest.split(','),
                                  args.pop_of_interest,
                                  args.pops_reference.split(','))
    elif args.type == 'sweep':
        merge_seeds_over_stp(args.read_path,
                             args.write_path,
                             args.selection_strengths.split(','),
                             [int(t) for t in args.times_of_sweep.split(',')],
                             args.pop_of_sweep,
                             args.pop_of_interest,
                             #args.pops_of_interest.split(','),
                             args.pops_reference.split(','),
                             sweep_pos=args.sweep_pos)


if __name__ == '__main__':
    main()
    # Run with python merge_compstats.py --read_path ../simulations/compstats/sweep/
    # --write_path ../training/simulations/sweep/ --pop_of_interest p1 --pops_reference p3 --type sweep
    # --times_of_sweep 10,15,20,25,30 --selection_strengths 0.015,0.025,0.035,0.045,0.06,0.08,0.1 --pop_of_sweep pa
