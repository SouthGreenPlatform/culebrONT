#!/bin/bash
#SBATCH --job-name CulebrONT
#SBATCH --output %x_%j.log
#SBATCH --error %x_%j.log

# adapt module load depending of your cluster
module load singularity/3.6.0
module load python/3.7.2
module load graphviz/2.40.1

# path to repertory containing config.yaml and cluster_config.yaml
CONFIG_DIR="/dir/to/CONFIG/"
# path to repertory where CulebrONT is intalled
CULEBRONT="/path/to/CulebrONT_pipeline/"

# SLURM JOBS
 snakemake  --nolock --use-conda --use-singularity --cores -p \
 -s $CULEBRONT/Snakefile  --latency-wait 600000 --keep-going --restart-times 0 --rerun-incomplete \
 --configfile $CONFIG_DIR/config.yaml \
  --cluster "python3 $CULEBRONT/slurm_wrapper.py $CONFIG_DIR/config.yaml $CONFIG_DIR/cluster_config.yaml" \
  --cluster-config $CONFIG_DIR/cluster_config.yaml \
  --cluster-status "python3 $CULEBRONT/slurm_status.py" \
  --conda-prefix $CULEBRONT/build_conda_envs
