#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --partition=long

export CULEBRONT="/home/orjuela/Documents/2021/CULEBRONT/MERGE/merge_v1.5.0/bkup/CulebrONT_pipeline"
profile=$CULEBRONT"/gift_files/CulebrONT_SLURM/"

submit_culebront.sh -c config.yaml -k cluster_config.yaml -p $profile