#!/bin/bash -l
#SBATCH --job-name=assembly_90_69_real
#SBATCH --output=ipyrad_output_%j
#SBATCH --error=ipyrad_err_%j
#SBATCH --time=160:00:00
#SBATCH --partition=fat
#SBATCH --clusters=cruncher
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200000
#SBATCH --exclusive

export PATH=/opt/miniconda3/bin/:$PATH
#export PATH=/data/home/wolfproj/wolfproj-03/miniconda3/bin:$PATH
ipyrad -p params-clean_90_69_real.txt -s 234567 -t 1 -c 20 -f