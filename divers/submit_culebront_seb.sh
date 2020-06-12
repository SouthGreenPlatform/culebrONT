#!/bin/bash
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --partition normal
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12

module load system/Miniconda3/1.0
#env=/home/$(whoami)/.conda/envs/snakemake
#[ ! -d $env ] && echo -e "## [$(date) - culebrONT]\t Creating conda environment for snakemake" && conda env create -f environment.yaml -n snakemake

source activate snakemake
module load system/singularity/3.3.0

configFile="config-hub_seb.yaml"


echo -e "## [$(date) - culebrONT]\t Starting main pipeline..."

snakemake --nolock --use-singularity --use-conda \
-p --verbose \
--latency-wait 60 --keep-going --restart-times 2 --rerun-incomplete \
--configfile ${configFile} \
--cores $SLURM_CPUS_PER_TASK \

# variales will be pass in arguments
#cluster_config="./cluster_config.yaml"
#datas_config="./config.yaml"
#scratch_dir="/scratch/"

#snakemake --nolock --use-singularity --use-conda --cores -p --verbose --keep-going --restart-times 2 --rerun-incomplete --latency-wait 60 -s Snakefile --jobs 4 --cluster "sbatch {cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu}{threads} " --cluster-config ${cluster_config} --configfile ${datas_config} --cluster-status "python3 slurm_status.py"

echo -e "## [$(date) - culebrONT]\t Creating snakemake report..."
snakemake   --configfile ${configFile} --report snakemake_report.html
