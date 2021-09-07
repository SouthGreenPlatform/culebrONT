#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --partition=long

export CULEBRONT="/path/to/CulebrONT_pipeline"
profile=$CULEBRONT"/gift_files/CulebrONT_SLURM/"

./submit_culebront.sh -c config.yaml -k cluster_config.yaml -p $profile