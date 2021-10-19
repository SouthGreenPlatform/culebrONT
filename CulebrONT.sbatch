#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --partition=long

#export culebront path to env
export CULEBRONT="/path/to/CulebrONT_pipeline"
#instead of this export, you can load the culebront module
#module load CulebrONT_pipeline-1.6.0/

profile=$CULEBRONT"/profiles/CulebrONT/"
#### Run a dry-run to verify install
# You can copy config.yaml and cluster_config.yaml from $CULEBRONT directly
# cp $CULEBRONT/config.yaml /path/to/
# cp $CULEBRONT/cluster_config.yaml /path/to/
submit_culebront.sh -c /path/to/config.yaml -k /path/to/cluster_config.yaml" -p $profile  -a  "--dryrun"