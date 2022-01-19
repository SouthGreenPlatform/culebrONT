#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --partition=long

module load culebrONT

culebrONT run_cluster -c config.yaml --dryrun