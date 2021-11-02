#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log
#SBATCH --partition=long

# CulebrONT needs 2 global variables : CULEBRONT and PROFILE
# if a module load was created, these variables were declared before, so load this module
module load CulebrONT/1.6.1

# if you don't use module load please export these two variables like this:
#export CULEBRONT="/path/to/CulebrONT_pipeline"
#export PROFILE="/path/to/your/own/profile/"

#### Run a dry-run to verify install
# You can copy config.yaml and cluster_config.yaml from $CULEBRONT directly and edit it
# cp $CULEBRONT/config.yaml /path/to/
# cp $CULEBRONT/cluster_config.yaml /path/to/

# Example in verbose mode and dry-run with default cluster_config included on profile
submit_culebront.sh -v -c /path/to/config.yaml -p $PROFILE  -a  "--dryrun"

# Example in verbose mode and dry-run with own cluster_config
submit_culebront.sh -v -c /path/to/config.yaml -k /path/to/cluster_config.yaml -p $PROFILE  -a  "--dryrun"

# Example of build dag or rulegraph (must disable verbose mode)
submit_culebront.sh -c /path/to/config.yaml -p $PROFILE  -a  "--dag" > schema_pipeline_samples.dot && dot schema_pipeline_samples.dot -Tpdf > schema_pipeline_samples.pdf
submit_culebront.sh -c /path/to/config.yaml -p $PROFILE  -a  "--rulegraph" > schema_pipeline_global.dot && dot schema_pipeline_global.dot -Tpdf > schema_pipeline_global.pdf

