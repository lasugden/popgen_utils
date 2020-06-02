import argparse
import pandas as pd
import os


def get_seeds(path):
    """
    From a directory, get all seed values by searching over names

    Args:
        path (str): path to search

    Returns:
        list of ints: the seeds that can be iterated over

    """
    seeds = os.listdir(path)
    seeds = list(set([int(x[:x.find('_')]) for x in seeds if x.find('100bins.norm') > -1]))
    return seeds


def read_files(directory_path, seed, pop_of_interest, pop_reference):
    """
    Read in IHS, XPEHH, ISAFE, and FST files, merge on position, and return all statistics

    Args:
        directory_path (str): path to directory containing files
        seed (int): seed used to create the files
        pop_of_interest (str): in the format 'pN' where N is an integer in 1-3
        pop_reference (str): in the format 'pN' where N is an integer in 1-3

    Returns:
        Pandas DataFrame: the merged output of the statistics

    """
    ihs_path = os.path.join(directory_path, '%i_%s.ihs.out.100bins.norm' % (seed, pop_of_interest))
    xpehh_path = os.path.join(directory_path, '%i_%s%s_W.xpehh.out' % (seed, pop_of_interest, pop_reference))
    isafe_path = os.path.join(directory_path, '%i_%s.iSAFE.out' % (seed, pop_of_interest))
    fst_path = os.path.join(directory_path, '%i_%s%s.weir.fst' % (seed, pop_of_interest, pop_reference))

    ihs_col, xpehh_col, isafe_col, fst_col = ('iHS', 'XP-EHH', 'Fst', 'iSAFE')

    if (os.path.exists(ihs_path) and
            os.path.exists(xpehh_path) and
            os.path.exists(isafe_path) and
            os.path.exists(fst_path)):
        # Read in files
        df_ihs = pd.read_csv(ihs_path, skiprows=0, header=None, delim_whitespace=True, usecols=range(7),
                             names=['locus_name', 'pos', 'freq', 'ihh1', 'ihh2', 'ihs', ihs_col])
        df_xpehh = pd.read_csv(xpehh_path, skiprows=1, header=None, delim_whitespace=True,
                               names=['locus_name', 'pos', 'gpos', 'p1', 'ihh1', 'p2', 'ihh2', xpehh_col])
        df_isafe = pd.read_csv(isafe_path, skiprows=1, header=None, delim_whitespace=True,
                               names=['pos', isafe_col, 'daf'])
        df_fst = pd.read_csv(fst_path, skiprows=1, header=None, delim_whitespace=True,
                             names=['chrom', 'pos', fst_col])

        # Match locations for 0-indexing and 1-indexing
        df_isafe['pos'] -= 1
        df_fst['pos'] -= 1

        # Merge dataframes
        df = df_ihs[['locus_name', 'pos', ihs_col]]\
            .merge(df_xpehh[['pos', xpehh_col]], how='outer', on='pos')\
            .merge(df_isafe[['pos', isafe_col]], how='outer', on='pos')\
            .merge(df_fst[['pos', fst_col]], how='outer', on='pos')

        # Fill in NaNs appropriately
        nan_names = df['locus_name'].isnull()
        df.loc[nan_names, 'locus_name'] = 'SNP_POS_' + df.loc[nan_names, 'pos'].apply(str)
        df = df.loc[df.isnull().sum(axis=1) < 4, :]
        df = df.fillna(-998)

        return df
    else:
        return None


def write_output(directory_path, seed, population, df):
    """
    Write the merged output to a file

    Args:
        directory_path (str): directory into which the file should be saved
        seed (int): the seed used to create the file type OR combination of metadata
        population (str): of the type pN where N is 1-2
        df (Pandas DataFrame): the combined statistics for writing

    """
    seed = '' if seed is None else seed
    file_name = '%s_%s_allstats.txt' % (str(seed), population)
    file_path = os.path.join(directory_path, file_name)
    df.to_csv(file_path, sep='\t', index=False)


def merge_all_seeds_and_write(input_directory,
                              output_directory,
                              pops_of_interest,
                              pop_reference):
    """
    Iterate over all seeds in a directory and merge.

    Args:
        input_directory (str): the directory in which all of the input data lives
        output_directory (str): the directory to which the combined statistics should be written
        pops_of_interest (str/list of str): in the format 'pN' where N is an integer in 1-3
        pop_reference (str): in the format 'pN' where N is an integer in 1-3

    """
    pops_of_interest = [pops_of_interest] if isinstance(pops_of_interest, str) else pops_of_interest
    seeds = get_seeds(input_directory)
    for seed in seeds:
        for population in pops_of_interest:
            df = read_files(input_directory, seed, population, pop_reference)
            if df is None:
                print('WARNING: some files are missing for seed %i and population %s' % (seed, population))
            else:
                write_output(output_directory, seed, population, df)


def merge_all_seeds_and_extract(input_directory,
                                output_directory,
                                pops_of_interest,
                                pop_reference,
                                sweep_pos,
                                additional_name_text=None):
    """
    Iterate over all seeds in a directory and merge.

    Args:
        input_directory (str): the directory in which all of the input data lives
        output_directory (str): the directory to which the combined statistics should be written
        pops_of_interest (str/list of str): in the format 'pN' where N is an integer in 1-3
        pop_reference (str): in the format 'pN' where N is an integer in 1-3
        sweep_pos (int): the position of the sweep to extract
        additional_name_text (str, optional): if set, this will be appended to filename

    """
    pops_of_interest = [pops_of_interest] if isinstance(pops_of_interest, str) else pops_of_interest
    seeds = get_seeds(input_directory)
    combined = []
    for population in pops_of_interest:
        for seed in seeds:
            df = read_files(input_directory, seed, population, pop_reference)
            if df is None:
                print('WARNING: some files are missing for seed %i and population %s' % (seed, population))
            else:
                combined.append(df.loc[df['pos'] == sweep_pos, :])

        if len(combined) > 0:
            df = pd.concat(combined, axis=0, ignore_index=True)
            write_output(output_directory, additional_name_text, population, df)


def merge_seeds_over_stp(input_directory,
                         output_directory,
                         selection_strength,
                         time_of_sweep,
                         pop_of_sweep,
                         pops_of_interest,
                         pop_reference,
                         sweep_pos=None):
    """
    Iterate over populations, selection strengths, and sweep times and extract combined statistics.

    Args:
        input_directory (str): OUTER path that includes all simulations
        output_directory (str): path into which allstats will be written
        selection_strength (float/list of floats): selection strength
        time_of_sweep (int/list of ints): kiloyears of sweep
        pops_of_interest (str/list of str): in the format 'pN' where N is an integer in 1-3
        pop_of_sweep (str): the population in which the sweep occured
        pop_reference (str): in the format 'pN' where N is an integer in 1-3
        sweep_pos (int, optional): the position of the sweep to extract
            if None, do not extract a single sweep_pos, but rather save each file

    """
    selection_strength = [selection_strength] if not isinstance(selection_strength, list) else selection_strength
    time_of_sweep = [time_of_sweep] if not isinstance(time_of_sweep, list) else time_of_sweep
    pops_of_interest = [pops_of_interest] if isinstance(pops_of_interest, str) else pops_of_interest

    for sel in selection_strength:
        for t in time_of_sweep:
            sub_input = os.path.join(input_directory, pop_of_sweep, 'T%i' % t, 's%s' % sel)
            if sweep_pos is None:
                merge_all_seeds_and_write(sub_input, output_directory, pops_of_interest, pop_reference)
            else:
                merge_all_seeds_and_extract(sub_input,
                                            output_directory,
                                            pops_of_interest,
                                            pop_reference,
                                            sweep_pos,
                                            't%i_s%s' % (t, sel))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_path', help='Path to the directory to iterate over', type=str)
    parser.add_argument('--write_path', help='Path to the directory to write into. Must already exist.', type=str)
    parser.add_argument('--pops_of_interest', help='Comma-separated list of pops (pN; with no spaces).', type=str)
    parser.add_argument('--pop_reference', help='Reference population (pN).', type=str)
    parser.add_argument('--type', help='Can be neutral or sweep', type=str, choices=['neutral', 'sweep'])
    parser.add_argument('--times_of_sweep', help='Comma-separated list of integer times of sweep', type=str)
    parser.add_argument('--selection_strengths', help='Comma-separated list of float selection strengths', type=str)
    parser.add_argument('--pop_of_sweep', help='Population of sweep (pN; with no spaces)', type=str)
    parser.add_argument('--sweep_pos', help='Sweep position', default=500000, type=int)

    args = parser.parse_args()

    if args.type == 'neutral':
        merge_all_seeds_and_write(args.read_path,
                                  args.write_path,
                                  args.pops_of_interest.split(','),
                                  args.pop_reference)
    elif args.type == 'sweep':
        merge_seeds_over_stp(args.read_path,
                             args.write_path,
                             args.selection_strengths.split(','),
                             [int(t) for t in args.times_of_sweep.split(',')],
                             args.pop_of_sweep,
                             args.pops_of_interest.split(','),
                             args.pop_reference,
                             sweep_pos=args.sweep_pos)


if __name__ == '__main__':
    main()
    # Run with python merge_compstats.py --read_path ../simulations/compstats/sweep/
    # --write_path ../training/simulations/sweep/ --pops_of_interest p1,p2 --pop_reference p3 --type sweep
    # --times_of_sweep 10,15,20,25,30 --selection_strengths 0.015,0.025,0.035,0.045,0.06,0.08,0.1 --pop_of_sweep pa
