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
from popgen_utils.haplotype_simulation import slim_population_definitions
from popgen_utils.haplotype_simulation import oscar_scripts


def slim_definition_to_input_neutral(definition_name,
                                     project_name,
                                     model_name,
                                     sims=1000,
                                     data_path=None):

    """
    Read in slim definition file and replace with neutral counts

    Args:
        definition_name (str): name of the definition file, e.g. gravel_model
        project_name (str): name of the project directory into which data was
            saved, e.g. 'sweep_hmm'
        model_name (str):the name of the particular model,
            e.g. 'gravel_neutral'
    """

   txt = ilresources.open_text(slim_population_definitions, f'{definition_name}.txt').read()

    # Make paths (e.g. datapath will be /users/la7/data/lalpert/data/)
    if data_path is None:
        data_path = config.params()['paths']['data']

    base_path = opath.join(data_path, project_name)
    if not opath.exists(base_path):
        os.mkdir(base_path)
    slim_path = opath.join(base_path, 'slim') 
    if not opath.exists(slim_path):
        os.mkdir(slim_path)
    slim_model_path = opath.join(slim_path, model_name)
    if not opath.exists(slim_model_path):
        os.mkdir(slim_model_path)

    # Save the metadata
    with open(opath.join(slim_path, f'{model_name}.yaml'), 'wb') as fp:
        yaml.dump({
            'date': datetime.now().strftime('%y%m%d'),
            'definition_name': definition_name,
            'project_name': project_name,
            'model_name': model_name,
            'sims':str(sims),
            'template_hash': hashing.hash(txt),
        }, fp, encoding='utf-8')

    for sim in range(int(sims)):
        parameter_model_name = (f'{model_name}_sim-{sim}')
        formatted_txt = txt.format(**{
            'vcf_file_output': opath.join(slim_model_path, f'{parameter_model_name}.vcf'),
            'ms_file_output': opath.join(slim_model_path, f'{parameter_model_name}_ms.txt'),
        })

        fp = open(opath.join(slim_model_path, f'{parameter_model_name}.slim'), 'w')
        fp.write(formatted_txt)
        fp.close()


def slim_definition_to_input(definition_name,
                             selection_coefficient_min,
                             selection_coefficient_max,
                             sweep_population,
                             sweep_time,
                             project_name,
                             model_name,
                             data_path=None,
                             sims_per_sweeptime=1000):
    """
    Read in a slim definition file and replace with the correct values from
    input data.

    Args:
        definition_name (str): name of the definition file, e.g. gravel_model
        selection_coefficient_min (float): the min selection
            coefficients to model, e.g. 0.02
        sweep_population (str or list of str): the populations to compare,
            e.g. ['p1', 'p2', 'p3']
        sweep_time (int or list of ints): number of kiloyears of sweep starts,
            e.g. [5, 10, 15, 20]
        project_name (str): name of the project directory into which data will
            be saved, e.g. 'sweep_hmm'
        model_name (str): the name of the particular model,
            e.g. 'gravel_sweep'
        data_path (str, optional): path into which the output will be saved.
            Defaults to the directory from config
        sims_per_sweeptime (int, optional): number of simulations per population and
            sweep time 

    """
    # Format the appropriate model file
    # noinspection PyTypeChecker
    txt = ilresources.open_text(slim_population_definitions, f'{definition_name}.txt').read()

    # Make paths (e.g. datapath will be /users/la7/data/lalpert/data/)
    if data_path is None:
        data_path = config.params()['paths']['data']

    base_path = opath.join(data_path, project_name)
    if not opath.exists(base_path):
        os.mkdir(base_path)
    slim_path = opath.join(base_path, 'slim') 
    if not opath.exists(slim_path):
        os.mkdir(slim_path)
    slim_model_path = opath.join(slim_path, model_name)
    if not opath.exists(slim_model_path):
        os.mkdir(slim_model_path)

        
    # Change this to draw selection coefficients at random from [selection_coefficient_min, selection_coefficient_max]
    # Make sure inputs are iterable
    #selection_coefficient = [selection_coefficient] if isinstance(selection_coefficient, float) \
       # else selection_coefficient
    selection_coefficient = [round(random.uniform(selection_coefficient_min,selection_coefficient_max), 8) for x in range(sims_per_sweeptime)]
    sweep_population = [sweep_population] if isinstance(selection_coefficient, str) else sweep_population
    sweep_time = [sweep_time] if isinstance(sweep_time, int) else sweep_time

    # Save the metadata
    with open(opath.join(slim_path, f'{model_name}.yaml'), 'wb') as fp:
        yaml.dump({
            'date': datetime.now().strftime('%y%m%d'),
            'definition_name': definition_name,
            'selection_coefficient': selection_coefficient,
            'sweep_population': sweep_population,
            'sweep_time': sweep_time,
            'project_name': project_name,
            'model_name': model_name,
            'template_hash': hashing.hash(txt),
        }, fp, encoding='utf-8')

    # Iterate
    for coeff in selection_coefficient:
        for pop in sweep_population:
            for time in sweep_time:
                # Convert sweep time into generations
                sweep_start_generation = 58000 - int(time*1000.0/25)

                parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                        f'pop-{pop}_start-{time}')

                # Format the text appropriately
                formatted_txt = txt.format(**{
                    'selection_coefficient': coeff,
                    'sweep_population': pop,
                    'sweep_start_generation': sweep_start_generation,
                    'vcf_file_output': opath.join(slim_model_path, f'{parameter_model_name}.vcf'),
                    'ms_file_output': opath.join(slim_model_path, f'{parameter_model_name}_ms.txt'),
                })

                fp = open(opath.join(slim_model_path, f'{parameter_model_name}.slim'), 'w')
                fp.write(formatted_txt)
                fp.close()


def run_slim_neutral(definition_name, project_name, model_name, data_path=None):
    """
    Run the slim binary multiple times for neutral simulations

    Args:
        definition_name (str): name of the definition file, e.g. gravel_model
        project_name (str): name of the project directory into which data was
            saved, e.g. 'sweep_hmm'
        model_name (str):the name of the particular model,
            e.g. 'gravel_neutral'
    """
 

    txt = ilresources.open_text(oscar_scripts, 'run_slim.txt').read()

    # Make paths (e.g. datapath will be /users/la7/data/lalpert/data/)
    if data_path is None:
        data_path = config.params()['paths']['data']

    base_path = opath.join(data_path, project_name)
    if not opath.exists(base_path):
        os.mkdir(base_path)
    slim_path = opath.join(base_path, 'slim') 
    if not opath.exists(slim_path):
        os.mkdir(slim_path)
    slim_model_path = opath.join(slim_path, model_name)
    if not opath.exists(slim_model_path):
        os.mkdir(slim_model_path)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)
    sims = params['sims']

    for sim in range(int(sims)):
        parameter_model_name = (f'{model_name}_sim-{sim}')
        formatted_txt = txt.format(**{
                        'slim_file':opath.join(slim_model_path,f'{parameter_model_name}.slim'),
                        'parameter_model_name': parameter_model_name,
        })
        fp = open(opath.join(slim_model_path, f'{parameter_model_name}.sh'), 'w')
        fp.write(formatted_txt)
        fp.close()
        os.system('sbatch '+opath.join(slim_model_path,f'{parameter_model_name}.sh'))



def run_slim(project_name, model_name, data_path=None):
    """
    Run the slim binary. Use the same names as those used for definition
    creation.

    Args:
        project_name (str): name of the project directory into which data was
            saved, e.g. 'sweep_hmm'
        model_name (str): the name of the particular model,
            e.g. 'gravel_sweep'

    Returns:

    """
    txt = ilresources.open_text(oscar_scripts, 'run_slim.txt').read()



    # Make paths (e.g. datapath will be /users/la7/data/lalpert/data/)
    if data_path is None:
        data_path = config.params()['paths']['data']

    base_path = opath.join(data_path, project_name)
    if not opath.exists(base_path):
        os.mkdir(base_path)
    slim_path = opath.join(base_path, 'slim') 
    if not opath.exists(slim_path):
        os.mkdir(slim_path)
    slim_model_path = opath.join(slim_path, model_name)
    if not opath.exists(slim_model_path):
        os.mkdir(slim_model_path)
    
    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)
    scoeffs = params['selection_coefficient']
    times = params['sweep_time']
    pops_of_interest = params['sweep_population']

    for pop in pops_of_interest:
        for time in times:
            for coeff in scoeffs:
                parameter_model_name = (f'{model_name}_coeff-{coeff}_pop-{pop}_start-{time}')
                if not os.path.isfile(opath.join(slim_model_path, f'{parameter_model_name}.vcf')) or not os.path.isfile(opath.join(slim_model_path, f'{parameter_model_name}_ms.txt')):

                    formatted_txt = txt.format(**{
                        'slim_file': opath.join(slim_model_path,f'{parameter_model_name}.slim'),
                        'parameter_model_name': parameter_model_name,
                    })
                    fp = open(opath.join(slim_model_path, f'{parameter_model_name}.sh'), 'w')
                    fp.write(formatted_txt)
                    fp.close()
                    os.system('sbatch '+opath.join(slim_model_path,f'{parameter_model_name}.sh'))


if __name__ == '__main__':
    slim_definition_to_input('gravel_model_100kb',
                             selection_coefficient_min=0.02,
                             selection_coefficient_max=0.12,
                             sweep_population=['p1', 'p2', 'p3'],
                             sweep_time=[5,10,15,20],
                             project_name='sweep_hmm',
                             model_name='gravel_sweep_100kb_uniform_scoeff')

    run_slim(project_name='sweep_hmm',
             model_name='gravel_sweep_100kb_uniform_scoeff')
