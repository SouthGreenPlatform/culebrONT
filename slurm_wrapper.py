#!/usr/bin/env python3
# written by Bao Tram Vi, modified by julie Orjuela (IRD)
from snakemake.logging import logger

#logger.info("wrapper to launch CulebrONT on slurm HPC")

'''
wrap to pass cluster configuration parameters to --cluster snake option in culebrONT
==>  ONLY TO SLURM !!
'''

import os
import sys
from snakemake.utils import read_job_properties
from snakemake import load_configfile

jobscript = sys.argv[-1]
config = sys.argv[1] 
cluster_config = sys.argv[2]
#logger.info(f"INFO: {jobscript} {config} {cluster_config}")

#read_job_propeties def reads the job properties defined in a snakemake jobscript and return a dict containing information about the job
job_properties = read_job_properties(jobscript)
config_properties = load_configfile(config)

cluster_properties = load_configfile(cluster_config)

rule = job_properties['rule']
jobid = job_properties['jobid']
log = rule

#logger.info("INFO job properties:")
#logger.info(job_properties)
#logger.info("INFO cluster properties : ")
#logger.info(cluster_properties)

#"cluster": {"cpus-per-task": 4, "ntasks": 1, "mem-per-cpu": "2", "partition": "normal", "output": "logs/stdout/run_flye/fastq=5percentB1-1", "error": "logs/error/run_flye/fastq=5percentB1-1"}}

# recovery wildcards in variables fastq, assemblers, busco_step 
try:
    fastq = job_properties['wildcards']['fastq']
    log += '_{}'.format(fastq)
except (AttributeError, KeyError):
    pass

try: 
    assemblers = job_properties['wildcards']['assemblers']
    log += '_{}'.format(assemblers)
except (AttributeError, KeyError):
    pass

try: 
    busco_step = job_properties['wildcards']['busco_step']
    log += '_{}'.format(busco_step)
except (AttributeError, KeyError):
    pass

# recovery cpu per task from dict properties
outdir = config_properties['DATA']['OUTPUT']
logdir = os.path.join(outdir, "slurm_log")
os.makedirs(logdir, exist_ok=True)

# def get_ressources(rule):
#     """
#     define ressources from cluster_config file or get params default for rule
#     """
#     if rule in cluster_properties and 'partition' in cluster_properties[rule]:
#        queue = cluster_properties[rule]['partition']
#        mempercpu = cluster_properties[rule]['mem-per-cpu']
#        cpus = cluster_properties[rule]['cpus-per-task']
#        return f"--partition {queue} --mem-per-cpu {mempercpu} --cpus-per-task {cpus} "
#     else:
#        queue = cluster_properties['__default__']['partition']
#        mempercpu = cluster_properties['__default__']['mem-per-cpu']
#        cpus = cluster_properties['__default__']['cpus-per-task']
#        return f'--partition {queue} --mem-per-cpu {mempercpu} --cpus-per-task {cpus} '
#

## getting ressources from cluster config given by user
def get_ressources(rule):
    if rule in cluster_properties and 'partition' in cluster_properties[rule]:
       queue = cluster_properties[rule]['partition']
       cpus = cluster_properties[rule]['cpus-per-task']
       if 'mem-per-cpu' in cluster_properties[rule].keys():
            mem = f"--mem-per-cpu {cluster_properties[rule]['mem-per-cpu']}"
       elif 'mem' in cluster_properties[rule].keys():
            mem = f"--mem {cluster_properties[rule]['mem']}"
       return f"--partition {queue} {mem} --cpus-per-task {cpus}"
    else:
       queue = cluster_properties['__default__']['partition']
       cpus = cluster_properties['__default__']['cpus-per-task']
       #mem = cluster_properties['__default__']['mem-per-cpu']
       if 'mem-per-cpu' in cluster_properties['__default__'].keys():
            mem = f"--mem-per-cpu {cluster_properties['__default__']['mem-per-cpu']}"
       elif 'mem' in cluster_properties['__default__'].keys():
            mem = f"--mem {cluster_properties['__default__']['mem']}"
       return f'--partition {queue} {mem} --cpus-per-task {cpus} '


partition = get_ressources(rule)

sbatch = f'sbatch --parsable --job-name {rule} {partition} --ntasks 1 --output {logdir}/{log}.log_%j --error {logdir}/{log}.log_%j'

# read jobscript and insert information such as node, used and memory
with open(jobscript, "r") as j:
    scripts = j.readlines()

scripts.insert(1, "echo -e \"# sbatch parameters: \"{}\"\"\n".format(sbatch))
scripts.insert(2, "echo -e \"# Job running on node: $SLURM_JOB_NODELIST\"\n")
scripts.insert(3, "echo -e \"# Number of used CPUS: $SLURM_CPUS_ON_NODE\"\n")
scripts.insert(4, "echo -e \"# Memory per CPU in megabyte: $SLURM_MEM_PER_CPU\"\n")
scripts.insert(5, "echo -e \"# Partition: $SLURM_JOB_PARTITION\"\n")

with open(jobscript, "w") as j:
    j.writelines(scripts)

cmdline = " ".join([sbatch, jobscript])
logger.info(f'INFO : {cmdline}')

os.system(cmdline)

