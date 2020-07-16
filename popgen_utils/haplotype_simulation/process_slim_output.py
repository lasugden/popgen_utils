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



def ms_to_mapfile_hapfile(filename, num_individuals, data_path=None):
    file = open(filename+'_ms.txt','r')
    f = file.read()
    file.close()
    f = f.strip().splitlines()
    sites = f[2].strip().split()[1:]
    sites = [float(x) for x in sites]

    start = 3
    pop_haplotypes = []
    for pop in range(len(num_individuals)):
        pop_haplotypes.append(f[start:start+num_individuals[pop]*2])
        start = start+num_individuals[pop]*2

    #p1ind = f[3:3+p1size*2]
    #p2ind = f[3+p1size*2:3+p1size*2+p2size*2]
    #p3ind = f[3+p1size*2+p2size*2:]
    out = open(filename+'_map.txt','w')
    SNPnum = 1
    outtext = ''
    oldmappos = 0
    oldphyspos = 0
    for mappos in sites:
        #mappos = site
        if mappos <= oldmappos:
            mappos = oldmappos + 0.000001
        oldmappos = mappos
        physpos = int(mappos*1e6)
        if physpos <= oldphyspos:
            physpos += 1
        oldphyspos = physpos
        outtext += '1 SNP'+str(SNPnum)+' '+str(mappos)+' '+str(physpos)+'\n'
        SNPnum += 1
    out.write(outtext.strip())
    out.close()

    for pop in range(len(num_individuals)):
        popname = 'p'+str(pop+1)
        out = open(filename+'_'+popname+'.hap','w')
        outtext = ''
        for line in pop_haplotypes[pop]:
            outtext += " ".join(line)+'\n'
        out.write(outtext.strip())
        out.close()


    # out = open(os.path.join(path,str(seed)+'_p1.hap'),'w')
    # outtext = ''
    # for line in p1ind:
    #     outtext += " ".join(line)+'\n'
    # out.write(outtext.strip())
    # out.close()

    # out = open(os.path.join(path,str(seed)+'_p2.hap'),'w')
    # outtext = ''
    # for line in p2ind:
    #     outtext += " ".join(line)+'\n'
    # out.write(outtext.strip())
    # out.close()

    # out = open(os.path.join(path,str(seed)+'_p3.hap'),'w')
    # outtext = ''
    # for line in p3ind:
    #     outtext += " ".join(line)+'\n'
    # out.write(outtext.strip())
    # out.close()


def slim_out_to_selscan(project_name, model_name, type, data_path=None):

    #number of individuals for each population -- match slim input files, made to match 1000G phase1 samples.
    #the number of haplotypes will be 2x the number of individuals
    #p1indiv = 108
    #p2indiv = 99
    #p3indiv = 103

    #numpops = 3
    num_individuals = [108, 99, 103]

    if data_path is None:
        data_path = config.params()['paths']['data']
    base_path = opath.join(data_path, project_name)
    slim_path = opath.join(base_path, 'slim')
    slim_model_path = opath.join(slim_path, model_name)

    yaml_file = open(opath.join(slim_path,f'{model_name}.yaml'))
    params = yaml.load(yaml_file)

    if type == 'sweep':
        scoeffs = params['selection_coefficient']
        times = params['sweep_time']
        pops_of_interest = params['sweep_population']
        for coeff in scoeffs:
            for time in times:
                for pop in pops_of_interest:
                    parameter_model_name = (f'{model_name}_coeff-{coeff}_'
                                        f'pop-{pop}_start-{time}')
                    filename = opath.join(slim_model_path, f'{parameter_model_name}')
                    if not opath.isfile(filename+'_p3.hap'):
                        ms_to_mapfile_hapfile(filename,num_individuals)

    if type == 'neutral':
        sims = params['sims']
        for sim in range(int(sims)):
            parameter_model_name = (f'{model_name}_sim-{sim}')
            filename = opath.join(slim_model_path,f'{oarameter_model_name}')
            if not opath.isfile(filename+'_p3.hap'):
                ms_to_mapfile_hapfile(filename,num_individuals)








