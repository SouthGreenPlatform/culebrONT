#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log

#module load system/Miniconda3/1.0

#env=/home/$(whoami)/.conda/envs/snakemake
#[ ! -d $env ] && echo -e "## [$(date) - culebrONT]\t Creating conda environment for snakemake" && conda env create -f environment.yaml -n snakemake

#source activate snakemake

module load system/singularity/3.3.0
module load system/python/3.7.2

snakemake --unlock

#snakemake --nolock --use-singularity --cores -p --verbose -s Snakefile-scp \
#--latency-wait 60 --keep-going --restart-times 2 --rerun-incomplete  \
#--cluster "python3 slurm_wrapper.py" \
#--cluster-config cluster_config.yaml \
#--cluster-status "python3 slurm_status.py"


snakemake --nolock --use-singularity --use-conda --cores -p --verbose -s Snakefile \
--latency-wait 60 --keep-going --restart-times 2 --rerun-incomplete \
--profile slurm-culebrONT \

#snakemake --nolock --use-singularity --use-conda --cores -p --verbose --keep-going --restart-times 2 --rerun-incomplete --latency-wait 60 -s Snakefile --jobs 100 --cluster "sbatch {cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu}{threads} " --cluster-config ${cluster_config} --configfile ${datas_config} --cluster-status "python3 slurm_status.py"

