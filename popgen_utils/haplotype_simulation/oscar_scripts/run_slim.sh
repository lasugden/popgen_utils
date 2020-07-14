#!/bin/bash
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mem=5G
#SBATCH -o '{parameter_model_name}.out'

module load slim
slim {slim_file}