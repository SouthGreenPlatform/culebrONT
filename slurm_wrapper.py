#!/usr/bin/env python3
# written by Bao Tram Vi, modified by julie Orjuela (IRD)

'''
wrap to pass cluster configuration parameters to --cluster snake option in culebrONT
'''
import os
import sys
from snakemake.utils import read_job_properties
from snakemake import load_configfile

jobscript = sys.argv[-1]
config = 'config.yaml'

#read_job_propeties def reads the job properties defined in a snakemake jobscript and return a dict containing information about the job
job_properties = read_job_properties(jobscript)
config_properties = load_configfile(config)
print(job_properties)


rule = job_properties['rule']
# jobid = job_properties['jobid']
log = rule

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
cpus_per_task = job_properties['threads']
# if int(cpus_per_task) < 8:
#     cpus_per_task = '8'
    
outdir = config_properties['DATA']['OUTPUT']
logdir = os.path.join(outdir, "slurm_log")
os.makedirs(logdir, exist_ok=True)

# ajouter ici la lecture du fichier cluster_config.yaml donné par l'utilisateur? 
#if rule == 'run_medaka_train':
#    partition = '--partition supermem --mem-per-cpu 4G'
#    cpus_per_task = '8'
#else:
#    partition = '--partition normal'

sbatch = f'sbatch --parsable --job-name {rule} {partition} --cpus-per-task {cpus_per_task} --ntasks 1 --output {logdir}/{log}.log_%j --error {logdir}/{log}.log_%j'

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
# sbatch --job-name {cluster.job-name} --partition {cluster.partition} --account {cluster.account} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output} --error {cluster.error}

os.system(cmdline)
print(cmdline)

