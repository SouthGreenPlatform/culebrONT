#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log

module load system/Miniconda3/1.0
env=${HOME}/.conda/envs/snakemake
[ ! -d $env ] && echo -e "## [$(date) - culebrONT]\t Creating conda environment for snakemake" && conda env create -f envs/environment.yaml -n snakemake

source activate snakemake
module load system/singularity/3.3.0
module load system/python/3.7.2

snakemake --unlock

# SLURM JOBS WITHOUT PROFILES
snakemake --nolock --use-conda --use-singularity --cores -p --verbose -s Snakefile \
--latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete  \
--configfile config.yaml \
--cluster "python3 slurm_wrapper.py config.yaml cluster_config.yaml" \
--cluster-config cluster_config.yaml \
--cluster-status "python3 slurm_status.py"

# USING PROFILES
#snakemake --nolock --use-singularity --use-conda --cores -p --verbose -s Snakefile --configfile config.yaml \
#--latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete --cluster-config cluster_config.yaml --profile slurm-culebrONT 


