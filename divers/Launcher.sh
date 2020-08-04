#!/bin/bash
#rm snakejob.*

# charging modules
module purge
module load system/python/3.7.2
module load bioinfo/Flye/2.4.1 
module load bioinfo/racon/1.4.3
module load bioinfo/minimap2/2.16

# Export enviromental variables
# DRMAA_LIBRARY_PATH autodetect wath scheduler is used by cluster
# to SGE only
# export DRMAA_LIBRARY_PATH=$SGE_ROOT/lib/lx-amd64/libdrmaa.so

# TODO : NT whereis DRMAA to SLURM variable $SLURM_ROOT?
# TODO:  --shadow-prefix ${scratch_dir} 

# building a graph to global pipeline 
snakemake -s Snakefile --configfile config.yaml --rulegraph | dot -Tpdf > schema_pipeline_global.pdf
# building a graph to pipeline by sample
snakemake -s Snakefile --configfile config.yaml --dag | dot -Tpdf > schema_pipeline_by_sample.pdf
# run with singularity
snakemake --use-singularity --latency 5555555 --rerun-incomplete

# variales will be pass in arguments 
cluster_config="./cluster_config.yaml"
datas_config="./config.yaml"
scratch_dir="/scratch/"

## pour SLURM config IRD
#snakemake  --latency-wait 5184000 -s Snakefile --jobs 100 --cluster "sbatch {cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu}{threads} " --cluster-config ${cluster_config} --configfile ${datas_config}

# pour SGE config IRD
snakemake  --latency-wait 5184000 -s Snakefile --jobs 100 --cluster "qsub {cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu}{threads} " --cluster-config ${cluster_config} --configfile ${datas_config}

#snakemake -s sRNA_pipeline.snake --filegraph | dot -Tpdf > schema_pipeline_files.pdf
#snakemake -s sRNA_pipeline.snake --dag | dot -Tpdf > schema_pipeline_samples.pdf
#snakemake -s sRNA_pipeline.snake --report REPORT.html

# add to clean files
#snakemake -s sRNA_pipeline.snake clean
