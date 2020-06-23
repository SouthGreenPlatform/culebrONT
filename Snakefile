#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import os, pprint
import re
from os import listdir
from snakemake.utils import validate
import warnings
from snakemake.logging import logger
from snakemake.dag import DAG

logger.info("""
    Welcome to CulebrONT  !
    Created on November 2019
    version: 1.0 July 2020
    @author: Julie Orjuela (IRD), Aurore Comte (IRD), Sebastien Ravel (CIRAD), Florian Charriat (CIRAD), Bao Tram Vi (IRD), François Sabot (IRD) and Sebastien Cunnac (IRD)
    @email: julie.orjuela@ird.fr, aurore@comte.ird.fr

#           ____   ____     _______     __________ ________   ________       ____    ___      _____________
#          6MMMMb/ `MM'     `M'`MM'     `MMMMMMMMM `MMMMMMMb. `MMMMMMMb.    6MMMMb   `MM\     `M'MMMMMMMMMM
#         8P    YM  MM       M  MM       MM      \  MM    `Mb  MM    `Mb   8P    Y8   MMM\     M /   MM   
#        6M      Y  MM       M  MM       MM         MM     MM  MM     MM  6M      Mb  M\MM\    M     MM
#        MM         MM       M  MM       MM    ,    MM    .M9  MM     MM  MM      MM  M \MM\   M     MM
#        MM         MM       M  MM       MMMMMMM    MMMMMMM(   MM    .M9  MM      MM  M  \MM\  M     MM
#        MM         MM       M  MM       MM    `    MM    `Mb  MMMMMMM9'  MM      MM  M   \MM\ M     MM
#        MM         MM       M  MM       MM         MM     MM  MM  \M\    MM      MM  M    \MM\M     MM
#        YM      6  YM       M  MM       MM         MM     MM  MM   \M\   YM      M9  M     \MMM     MM
#         8b    d9   8b     d8  MM    /  MM      /  MM    .M9  MM    \M\   8b    d8   M      \MM     MM
#          YMMMM9     YMMMMM9  _MMMMMMM _MMMMMMMMM _MMMMMMM9' _MM_    \M\_  YMMMM9   _M_      \M    _MM_


    Please cite our github https://github.com/SouthGreenPlatform/CulebrONT_pipeline
    Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) 
    and GPLv3 Intellectual property belongs to IRD and authors.
    """)

################ CONFIGFILE EXTRACTION ################
def install_dag_hook(callback):
    def postprocess_hook(self, __origmeth=DAG.postprocess):
        __origmeth(self)
        callback(self)
    DAG.postprocess = postprocess_hook

def dag_finalized(dag):
    for j in dag.jobs:
        print(j.output, j.output_mintime, j.input_maxtime)

### to use on V2
#install_dag_hook(dag_finalized)
#pprint.pprint (DAG.__dict__)

# importing configuration files
#pprint.pprint(workflow.__dict__)
if len(workflow.overwrite_configfiles) == 0:
    logger.info("You need to use --configfile option to snakemake command line")
    raise ValueError("You have to use --configfile option to snakemake command line")
else:
    path_config = workflow.overwrite_configfiles[0]

#configfile:'config.yaml'
cluster_config: 'cluster_config.yaml'

# --- Verification Configuration Files --- #
def get_fastq_list():
    list_of_fq = []
    for f in listdir(config['DATA']['FASTQ']):
        filename, file_extension = os.path.splitext(f)
        if file_extension in [".fastq",".fq",".fq.gz",".fastq.gz"]:
            list_of_fq.append(filename)
    return list_of_fq

# Configuration file checking
def check_config_dic(config):
    if re.match(r"\S+\/$", config["DATA"]["FASTQ"]) == None:
        config["DATA"]["FASTQ"] = f"{config['DATA']['FASTQ']}/"
    if not os.path.exists(config["DATA"]["FASTQ"]):
        logger.info("CONFIG FILE CHECKING FAIL : in the DATA section, FASTQ directory does not exists")
        raise ValueError("CONFIG FILE CHECKING FAIL : in the DATA section, FASTQ directory does not exists")
    if re.match(r"\S+\/$", config["DATA"]["FAST5"]) == None:
        config["DATA"]["FAST5"] = f"{config['DATA']['FAST5']}/"
    if not os.path.exists(config["DATA"]["FAST5"]):
        os.makedirs(config["DATA"]["FAST5"])
    else:
        fastq_list = get_fastq_list()
        fast5_list = os.listdir(config["DATA"]["FAST5"])
        diff = set(fastq_list) - set(fast5_list)
        if diff:
            warnings.warn(f"\n CONFIG FILE CHECKING WARNING : You don't have a fast5 repository for each of your fastq file (they should have the same name). This can raise a problem if you choose to use Nanopolish. Please check your data.")
            logger.info("CONFIG FILE CHECKING WARNING : You don't have a fast5 repository for each of your fastq file (they should have the same name). This can raise a problem if you choose to use Nanopolish. Please check your data.")
    if re.match(r"\S+\/$", config["DATA"]["ILLUMINA"]) == None:
        config["DATA"]["ILLUMINA"] = f"{config['DATA']['ILLUMINA']}/"
        if not os.path.exists(config["DATA"]["ILLUMINA"]):
            os.makedirs(config["DATA"]["ILLUMINA"])
    if re.match(r"\S+\/$", config["DATA"]["OUTPUT"]) == None:
        config["DATA"]["OUTPUT"] = f"{config['DATA']['OUTPUT']}/"
    if not os.path.isfile(config["DATA"]["REF"]):
        logger.info("CONFIG FILE CHECKING FAIL : in the DATA section, the REF file is not a file")
        raise ValueError("CONFIG FILE CHECKING FAIL : in the DATA section, the REF file is not a file")
    if config["QUALITY"]["ASSEMBLY"] == False and config["QUALITY"]["POLISHING"] == False and config["QUALITY"]["CORRECTION"] == False:
        logger.info("CONFIG FILE CHECKING FAIL : you need to set True for at least one quality step (ASSEMBLY, POLISHING or CORRECTION)")
        raise ValueError("CONFIG FILE CHECKING FAIL : you need to set True for at least one quality step (ASSEMBLY, POLISHING or CORRECTION)")
    else:
        if config["ASSEMBLY"]["CANU"] == False and config["ASSEMBLY"]["FLYE"] == False and config["ASSEMBLY"]["MINIASM"] == False:
            logger.info("CONFIG FILE CHECKING FAIL : you need to set True for at least one Assembly tool (CANU, FLYE, MINIASM ...")
            raise ValueError("CONFIG FILE CHECKING FAIL : you need to set True for at least one Assembly tool (CANU, FLYE, MINIASM ...)")

if (validate(config, "schemas/config.schema.yaml")) is not None:
    logger.info("CONFIG FILE CHECKING STRUCTURE FAIL : you need to verify confi.yaml KEYS:VALUES")
    raise ValueError("CONFIG FILE CHECKING STRUCTURE FAIL : you need to verify confi.yaml KEYS:VALUES.)")

check_config_dic(config)

# Cleaning Repport and dags if we rerun the snakemake a second time

def cleaning_for_rerun(config):
    if os.path.exists(f"{config['DATA']['OUTPUT']}/REPORT") or os.path.isfile(f"{config['DATA']['OUTPUT']}/dag.png"):
        logger.info("If you want to rerun culebront a second time please delete the {config['DATA']['OUTPUT']}REPORT directory and the {config['DATA']['OUTPUT']}dag.png file.")
        warnings.warn(f"\n If you want to rerun culebront a second time please delete the {config['DATA']['OUTPUT']}REPORT directory and the {config['DATA']['OUTPUT']}dag.png file. \n")

cleaning_for_rerun(config)

# Getting paths on usefully variables
output_dir = Path(config['DATA']['OUTPUT']).resolve().as_posix() + "/"
fastq = Path(config['DATA']['FASTQ']).resolve().as_posix()
ref = Path(config['DATA']['REF']).resolve().as_posix()
fast5 = Path(config['DATA']['FAST5']).resolve().as_posix()
illumina = Path(config['DATA']['ILLUMINA']).resolve().as_posix() + "/"
path_config = Path('config.yaml').resolve().as_posix()
path_snake = Path('Snakefile').resolve().as_posix()

#print (type(path_config))
#print (type(illumina))

def getting_ext(chaine):
    """
    create a Path object and return first suffix found from fist file
    """
    basepath = Path(chaine)
    for entry in basepath.iterdir():
        if entry.is_file():
            #print (entry.name)
            return entry.suffix

ext_illumina = getting_ext(illumina)
#print ('**********', ext_illumina)


# Declaring model variable to medaka
if 'model' in config['DATA'].keys():
    model = True
else:
    model = False

################ WILDCARDS  ################
regex=r'.*q\.{0,1}(gz){0,1}' # '.*(q|gz)$'
FASTQ, EXT, = glob_wildcards(f"{config['DATA']['FASTQ']}{{fastq}}.{{ext,{regex}}}")# TODO : mettre obligatoirement des fichiers en fastq.gz
nb = str(config['params']['RACON']['RACON_ROUNDS'])

# controling wildcard from config file params
ASSEMBLY_TOOLS = []
if config['ASSEMBLY']['CANU']:
    ASSEMBLY_TOOLS.append("CANU")
if config['ASSEMBLY']['FLYE']:
    ASSEMBLY_TOOLS.append("FLYE")
if config['ASSEMBLY']['MINIASM']:
    ASSEMBLY_TOOLS.append("MINIASM")
#print (ASSEMBLY_TOOLS)

CORRECTION_TOOLS = []
if config['CORRECTION']['MEDAKA']:
    CORRECTION_TOOLS.append("MEDAKA")
if config['CORRECTION']['NANOPOLISH']:
    CORRECTION_TOOLS.append("NANOPOLISH")
#print (CORRECTION_TOOLS)

POLISHING_TOOLS = ["RACON"]
# if config['POLISHING']['RACON']:
#     POLISHING_TOOLS.append("RACON")
#print (POLISHING_TOOLS)

BUSCO_STEPS = []
if config['QUALITY']['ASSEMBLY']:
    BUSCO_STEPS.append("STEP_ASSEMBLY")
if config['QUALITY']['POLISHING']:
    for polisher in POLISHING_TOOLS:
        BUSCO_STEPS.append("STEP_POLISHING_" + polisher)
if config['QUALITY']['CORRECTION']:
    for corrector in CORRECTION_TOOLS:
        BUSCO_STEPS.append("STEP_CORRECTION_" + corrector)

# Make sure running Mauve makes sense
if len(ASSEMBLY_TOOLS) <2 and len(BUSCO_STEPS) < 2 and config['MSA']['MAUVE']:
    warnings.warn("MAUVE is irrelevant if you have a single final assembly as your config file implies (need more than one assembler and/or a correction step). Mauve will not be run !! ")
    config['MSA']['MAUVE'] = False

# Make sure running fixstart makes sense
if not config['DATA']['CIRCULAR']:
    warnings.warn(f"\n Fixstart is irrelevant if you have not activated CIRCULAR. Fixstart will not be run !! \n")
    config['MSA']['FIXSTART'] = False

#at least a quality step must to be activated by user in  config
if len(BUSCO_STEPS) == 0 :
    raise ValueError("you have to activate at least a QUALITY (ASSEMBLY, POLISHING or CORRECTION) step in config file")

if len(BUSCO_STEPS) == 0 :
    raise ValueError("you have to activate at least a QUALITY (ASSEMBLY, POLISHING or CORRECTION) step in config file")

if int(config['params']['QUAST']['GENOME_SIZE_PB']) >= 100000000:
    warnings.warn(f"\n Weesam if fixed to FALSE because genome size !! \n")
    config['QUALITY']['WEESAM'] = False

if config['DATA']['CIRCULAR']:
    add_circular_name = "CIRCULARISED"
else:
    add_circular_name = ""

############################### DEF ####################################

draft_to_correction = f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{nb}/assembly.racon{nb}.fasta" if int(nb)>=1 else f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/Assembly{add_circular_name}.fasta"
draft_to_correction_index_fai = f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{nb}/assembly.racon{nb}.fasta.fai" if int(nb)>=1 else f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/Assembly{add_circular_name}.fasta.fai"
draft_to_correction_index_mmi = f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{nb}/assembly.racon{nb}.fasta.mmi" if int(nb)>=1 else f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/Assembly{add_circular_name}.fasta.mmi"

def get_threads(rule, default):
    """
    use threads define in cluster_config rule or rule default or default in snakefile
    """
    if rule in cluster_config and 'threads' in cluster_config[rule]:
        return int(cluster_config[rule]['threads'])
    elif '__default__' in cluster_config and 'threads' in cluster_config['__default__']:
        return int(cluster_config['__default__']['threads'])
    return default

# def get_fastq(wildcards):
#     for f in listdir(config['DATA']['FASTQ']):
#         filename, file_extension = os.path.splitext(f)
#         if filename == wildcards.fastq:
#             return f"{config['DATA']['FASTQ']}{f}"

def get_fastq(wildcards):
    for f in listdir(config['DATA']['FASTQ']):
        # filename, file_extension = os.path.splitext(f)
        filename = f.split('.')[0]
        if filename == wildcards.fastq:
            return f"{config['DATA']['FASTQ']}{f}"

def draft_to_racon(wildcards):
    n = int(wildcards.nb)
    if n == 1:
        if config['DATA']['CIRCULAR']:
            return f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/ROTATE/rotate_{n}/assembly.racon{n}.fasta"
        else:
            return f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/assembly{add_circular_name}.fasta"
    elif n > 1:
        if config['DATA']['CIRCULAR']:
            return f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/ROTATE/rotate_{n}/assembly.racon{n}.fasta"
        else:
            return f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{n-1}/assembly.racon{n-1}.fasta"
    else:
        raise ValueError(f"loop numbers must be 1 or greater: received {n}")

def draft_to_rotate(wildcards):
    n = int(wildcards.nb)
    if n == 1:
        return f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/assemblyCIRCULARISED_circTag.fasta"
    elif n > 1:
        #rules.run_racon.output.fasta
        return f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{n-1}/assembly.racon{n-1}.fasta"
    else:
        raise ValueError(f"loop numbers must be 1 or greater: received {n}")

def input_last():
    liste=[]
    #if len(BUSCO_STEPS) == 0 :
        #return liste.append("")
    if len(CORRECTION_TOOLS)>=1:
        for toolcorrection in CORRECTION_TOOLS:
            liste.append("STEP_CORRECTION_" + toolcorrection)
        return liste
    if len(POLISHING_TOOLS)>=1:
        for toolpolish in POLISHING_TOOLS:
            liste.append("STEP_POLISHING_" + toolpolish)
        return liste
    if len(ASSEMBLY_TOOLS)>=1:
        liste.append("STEP_ASSEMBLY")
        return liste
    else:
        raise ValueError("you have to activate at least a QUALITY (ASSEMBLY, POLISHING or CORRECTION) step in config file")

def output_final(wildcards):
    dico_final = {
        #"dag": f"{output_dir}dag.png"
    }
    if config['QUALITY']['BLOBTOOLS']:
        dico_final.update({
            "blob_files": expand(f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/BLOBTOOLS/output.quality.blobDB.table.txt", fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last())
        })
    if config['QUALITY']['WEESAM']:
        dico_final.update({
             "weesam_files": expand(f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/WEESAM/minimap2mapping.txt", fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last())
        })
    if config['QUALITY']['ASSEMBLYTICS']:
        dico_final.update({
             "assemblytics_files": expand(f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/ASSEMBLYTICS/OUT.Assemblytics_structural_variants.summary", fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last())
        })
    if config['QUALITY']['KAT']:
        dico_final.update({
            "kat_files":  expand(f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.gcp", fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last())
        })
    if config['DATA']['CIRCULAR'] and config['MSA']['FIXSTART'] and not config['MSA']['MAUVE'] :
        dico_final.update({
            "fixstart_files": expand(f"{output_dir}{{fastq}}/{{assemblers}}/MSA/FIXSTART-{{busco_step}}/startfixed_asm.fasta", fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last())
        })
    if config['MSA']['MAUVE']:
        dico_final.update({
            "mauve_files": expand(f"{output_dir}{{fastq}}/MAUVE_ALIGN/candidate_assemblies.xmfa", fastq=FASTQ)
        })
    return dico_final

def output_values_final(wildcards):
    dico = output_final(wildcards)
    val = list()
    for valeur in dico.values():
        val.append(str(valeur[0]))
    #print(val)
    return val

################################ FINAL ####################################
rule final:
    input:
        expand(f"{output_dir}REPORT/Report.html", fastq=FASTQ)

rule run_dico_final:
    input:
        unpack(output_final)
        #output_values_final
    output:
        check = f'{output_dir}ok.txt'
    log:
        error = f"{output_dir}LOGS/DICO-FINAL/dico-final.e"
    message:
        """
        check = {output.check}
        """
    shell:
        """
        echo 'Thanks for using CulebrONT' 1>{output.check} 2>{log.error}
        """
################################ ASSEMBLY ####################################
rule run_flye:
    """
    launch flye
    """
    threads: get_threads('run_flye', 8)
    input:
        fastq = get_fastq,
    output:
        fasta = f"{output_dir}{{fastq}}/FLYE/ASSEMBLER/assembly{add_circular_name}.fasta",
        info = f"{output_dir}{{fastq}}/FLYE/ASSEMBLER/assembly_info.txt",
    params:
        fasta_dir = directory(f"{output_dir}{{fastq}}/FLYE/ASSEMBLER/"),
        genome_size = config['DATA']['GENOME_SIZE'],
        circular = "--plasmids" if config['DATA']['CIRCULAR'] else "",
        move = "mv " if config['DATA']['CIRCULAR'] else "echo ",
    log:
        output = f"{output_dir}LOGS/ASSEMBLER/FLYE/{{fastq}}_FLYE.o",
        error = f"{output_dir}LOGS/ASSEMBLER/FLYE/{{fastq}}_FLYE.e",
    benchmark:
        f"{output_dir}LOGS/ASSEMBLER/FLYE/{{fastq}}_FLYE-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads: {threads}
        input:
            fastq: {input.fastq}
        output:
            fasta: {output.fasta}
            info: {output.info}
        params:
            fasta_dir: {params.fasta_dir}
            genome_size:   {params.genome_size}
            circular:  {params.circular}
        log:
            output: {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['FLYE_SIMG']
    shell:
        """
        flye --nano-raw {input.fastq} --genome-size {params.genome_size} {params.circular} --out-dir {params.fasta_dir} --threads {threads} 1>{log.output} 2>{log.error}
        {params.move} {params.fasta_dir}/assembly.fasta {output.fasta} 1>>{log.output} 2>>{log.error}
        """

rule run_canu:
    """
    launch canu
    """
    threads: get_threads('run_canu', 8)
    input:
        fastq = get_fastq,
    output:
        fasta_canu = f"{output_dir}{{fastq}}/CANU/ASSEMBLER/out_canu.contigs.fasta",
        fasta = f"{output_dir}{{fastq}}/CANU/ASSEMBLER/assembly2Circ.fasta" if config['DATA']['CIRCULAR'] else f"{output_dir}{{fastq}}/CANU/ASSEMBLER/assembly.fasta",
        trim_corr_fq = f"{output_dir}{{fastq}}/CANU/ASSEMBLER/out_canu.trimmedReads.fasta.gz"
    params:
        out_dir = directory(f"{output_dir}{{fastq}}/CANU/ASSEMBLER/"),
        genome_size = f"{config['DATA']['GENOME_SIZE']}",
        max_memory = f"{config['params']['CANU']['MAX_MEMORY']}",
        options = f"{config['params']['CANU']['OPTIONS']}",
    log:
        output = f"{output_dir}LOGS/ASSEMBLER/CANU/{{fastq}}_CANU.o",
        error = f"{output_dir}LOGS/ASSEMBLER/CANU/{{fastq}}_CANU.e",
    benchmark:
        f"{output_dir}LOGS/ASSEMBLER/CANU/{{fastq}}_CANU-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads: {threads}
        input:
            fastq: {input.fastq}
        output:
            fasta_canu: {output.fasta_canu}
            fasta: {output.fasta}
            trim_forr_fq: {output.trim_corr_fq}
        params:
            out_dir: {params.out_dir}
            genome_size: {params.genome_size}
            maxMemory: {params.max_memory}
            options: {params.options}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['CANU_SIMG']
    shell:
        """
        canu -p out_canu useGrid=false maxThreads={threads} maxMemory={params.max_memory} -d {params.out_dir} genomeSize={params.genome_size} -nanopore-raw {input.fastq} {params.options} 1>{log.output} 2>{log.error}
        ln -s {output.fasta_canu} {output.fasta} 1>>{log.output} 2>>{log.error}
        """

rule run_miniasm:
    """
    launch miniasm
    """
    threads: get_threads('run_miniasm', 4)
    input:
        fastq = get_fastq,
    output:
        gfa_miniasm = f"{output_dir}{{fastq}}/MINIASM/ASSEMBLER/output_miniasm.gfa", # voir si à utiliser par gfapy
    params:
        temp_paf = temp(f"{output_dir}{{fastq}}/MINIASM/ASSEMBLER/output_minimap2.paf"),
    log:
        error = f"{output_dir}LOGS/ASSEMBLER/MINIASM/{{fastq}}_MINIASM.e",
    benchmark:
        f"{output_dir}LOGS/ASSEMBLER/MINIASM/{{fastq}}_MINIASM-BENCHMARK.txt"
    priority: 30
    message:
           """
           Launching {rule}
           threads : {threads}
           input:
               fastq : {input.fastq}
           params:
               paf : {params.temp_paf}
           output:
               gfa : {output.gfa_miniasm}
           log:
               error: {log.error}
           """
    singularity:
        config['tools']['MINIASM_SIMG']
    shell:
         """
         minimap2 -x ava-ont -t {threads} {input.fastq} {input.fastq} 1> {params.temp_paf} 2>{log.error}
         miniasm -f {input.fastq} {params.temp_paf} 1> {output.gfa_miniasm} 2>>{log.error}
         """
    #[ -s {params.temp_paf} ] || echo "{params.temp_paf} is empty" 2>>{log.error}

rule run_minipolish:
    """
    launch minipolish
    """
    threads: get_threads('run_minipolish', 4)
    input:
        fastq = get_fastq,
        gfa_miniasm = rules.run_miniasm.output.gfa_miniasm
    output:
        gfa_minipolish = f"{output_dir}{{fastq}}/MINIASM/ASSEMBLER/output_minipolish.gfa",
        fasta = f"{output_dir}{{fastq}}/MINIASM/ASSEMBLER/assembly{add_circular_name}.fasta", #### {add_circular_name}
        info = f"{output_dir}{{fastq}}/MINIASM/ASSEMBLER/assembly_info.txt"
    params:
        #racon_rounds = config['params']['MINIPOLISH']['RACON_ROUNDS'],
        racon_rounds = "2",
    log:
        error = f"{output_dir}LOGS/ASSEMBLER/MINIASM/{{fastq}}_MINIASM_MINIPOLISH.e"
    benchmark:
        f"{output_dir}LOGS/ASSEMBLER/MINIASM/{{fastq}}_MINIASM_MINIPOLISH-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
            fastq : {input.fastq}
            gfa : {input.gfa_miniasm}
        output:
            gfa : {output.gfa_minipolish}
            fasta : {output.fasta}
            info : {output.info}
        params:
            racon_rounds: {params.racon_rounds}
        log:
            error: {log.error}
        """
    singularity:
        config['tools']['MINIPOLISH_SIMG']
    shell:
        """
        minipolish -t {threads} --rounds {params.racon_rounds} {input.fastq} {input.gfa_miniasm} 1> {output.gfa_minipolish} 2>{log.error}
        awk '/^S/{{print \">\"$2\"\\n\"$3}}' {output.gfa_minipolish} | fold > {output.fasta}
        ln -s {output.gfa_minipolish} {output.info} 2>>{log.error}
        """

############################### CIRCULARISATION POST CANU ##############################

rule run_circlator:
    """
    launch Circlator
    """
    threads: get_threads('run_circlator', 4)
    input:
        draft = rules.run_canu.output.fasta,
        fastq = rules.run_canu.output.trim_corr_fq,
        #fastq = get_fastq,
    output:
        fasta = f"{output_dir}{{fastq}}/CANU/ASSEMBLER/assemblyCIRCULARISED.fasta",
        info = f"{output_dir}{{fastq}}/CANU/ASSEMBLER/assembly_info.txt",
    params:
        log_mv = f"{output_dir}{{fastq}}/CANU/ASSEMBLER/circlator.log",
        out_dir = directory(f"{output_dir}{{fastq}}/CANU/ASSEMBLER/CIRCLATOR/"),
        options = f"{config['params']['CIRCLATOR']['OPTIONS']}",
    log:
        output = f"{output_dir}LOGS/ASSEMBLER/CIRCLATOR/{{fastq}}_CANU.o",
        error = f"{output_dir}LOGS/ASSEMBLER/CIRCLATOR/{{fastq}}_CANU.e",
    benchmark:
        f"{output_dir}LOGS/ASSEMBLER/CIRCLATOR/{{fastq}}_CANU-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
             fastq : {input.fastq}
             draft : {input.draft}
        output:
             fasta : {output.fasta}
             info : {output.info}
        params:
             out_dir: {params.out_dir}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        #config['tools']['CIRCLATOR_SIMG']
        config['tools']['MINICONDA_SIMG']
    conda: "envs/run_circlator_cenv.yml"
    shell:
        """
        rm -rf {params.out_dir} 1>{log.output} 2>{log.error}
        circlator all --thread {threads} {params.options} --verbose --bwa_opts "-x ont2d" {input.draft} {input.fastq} {params.out_dir}  1>>{log.output} 2>>{log.error}
        mv {params.out_dir}06.fixstart.fasta {output.fasta} 1>>{log.output} 2>>{log.error}
        mv {params.out_dir}06.fixstart.log {params.log_mv} 1>>{log.output} 2>>{log.error}
        ln -s "{params.out_dir}04.merge.circularise.log" "{output.info}" 1>>{log.output} 2>>{log.error}
        """

############################### TAGGING OF CIRCULAR MOLECULES ##########################

rule tag_circular:
    """
    Tagging title of circular molecules in assembly fasta files
    """
    threads: get_threads('tag_circular', 1)
    input:
        fasta = f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/assemblyCIRCULARISED.fasta",
        info = f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/assembly_info.txt",
    output:
        tagged_fasta = f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/assemblyCIRCULARISED_circTag.fasta"
    log:
        output = f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/tag_circular.log"
        #output = f"{output_dir}LOGS/ASSEMBLER/CIRCLATOR/{{fastq}}_CANU.o",
        #error = f"{output_dir}LOGS/ASSEMBLER/CIRCLATOR/{{fastq}}_CANU.e",
    benchmark:
        f"{output_dir}LOGS/TAGGING/{{fastq}}_{{assemblers}}-TAG-CIRCULAR-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
            fasta : {input.fasta}
            info : {input.info}
        output:
            tagged_fasta : {output.tagged_fasta}
        log:
            output = {log.output}",
        """
    singularity:
        config['tools']['MINICONDA_SIMG']
    conda: "envs/R_for_culebront_cenv.yml"
    shell:
        """
        exec > >(tee "{log.output}") 2>&1
        Rscript AdditionalScripts/tagCircSeq.R --seqFile="{input.fasta}" --logFile="{input.info}" --outFilePath="{output.tagged_fasta}"
        """

################################ INDEXING ###################################

rule index_fasta_to_correction:
    """
    create a .fai and a .mmi for each assembly fasta
    """
    threads: get_threads('index_fasta_to_correction', 4)
    input:
        draft = draft_to_correction
    output:
        index_fai = draft_to_correction_index_fai,
        index_mmi = draft_to_correction_index_mmi
    singularity:
        config['tools']['MEDAKA_SIMG']
    priority: 30
    shell:
        """
            samtools faidx {input.draft}
            minimap2 -d {input.draft}.mmi {input.draft}
        """

################################ POLISHING ####################################

rule run_racon:
    """
    launch Racon recursively n times (given by config.yaml)
    """
    threads: get_threads('run_racon', 4)
    input:
        draft = draft_to_racon,
        fastq = get_fastq,
    output:
        paf = f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{{nb}}/assembly.minimap4racon{{nb}}.paf",
        fasta = f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{{nb}}/assembly.racon{{nb}}.fasta"
    wildcard_constraints:
        nb = "[0-9]"
    log:
        output=f"{output_dir}LOGS/POLISHING/RACON/{{fastq}}_{{assemblers}}_RACON{{nb}}.o",
        error = f"{output_dir}LOGS/POLISHING/RACON/{{fastq}}_{{assemblers}}_RACON{{nb}}.e"
    benchmark:
        f"{output_dir}LOGS/POLISHING/RACON/{{fastq}}_{{assemblers}}_RACON{{nb}}-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
            draft : {input.draft}
            fastq : {input.fastq}
        output:
            paf : {output.paf}
            fasta : {output.fasta}
        log:
            error: {log.error}
        """
    singularity:
        config['tools']['RACON_SIMG']
    shell:
        """
        minimap2 -t {threads} {input.draft} {input.fastq} 1> {output.paf} 2>{log.error}
        racon -t {threads} {input.fastq} {output.paf} {input.draft} 1> {output.fasta} 2>>{log.error}
        """

rule rotate_circular:
    """
    Rotate circular
    """
    threads: get_threads('rotate_circular', 4)
    input:
        draft = draft_to_rotate,
    output:
        rotated = f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/ROTATE/rotate_{{nb}}/assembly.racon{{nb}}.fasta"
    log:
        output=f"{output_dir}LOGS/ROTATE/{{fastq}}_{{assemblers}}_RACON{{nb}}.o",
        error = f"{output_dir}LOGS/ROTATE/{{fastq}}_{{assemblers}}_RACON{{nb}}.e"
    benchmark:
        f"{output_dir}LOGS/ROTATE/{{fastq}}_{{assemblers}}_RACON{{nb}}-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
            draft : {input.draft}
        output:
            rotated : {output.rotated}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['MINICONDA_SIMG']
    conda:
        "envs/R_for_culebront_cenv.yml"
    shell:
        """
        Rscript AdditionalScripts/rotateCircSeqs.R --seqFile "{input.draft}" --outFilePath "{output.rotated}" 2>{log.error}
        """

rule run_nanopolish :
    """
    launch makerange, consensus and vcf2fasta
    """
    threads: get_threads('run_nanopolish', 8)
    input:
        draft = draft_to_correction,
    output:
        temp_bam = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/NANOPOLISH/reads.sorted.bam",
        temp_fasta = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/NANOPOLISH/reads.fasta",
        fasta = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/NANOPOLISH/consensus.fasta",
    params:
        fastq = get_fastq,
        fast5 = f"{config['DATA']['FAST5']}/{{fastq}}",
        directory = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/NANOPOLISH",
        segment = config['params']['NANOPOLISH']['NANOPOLISH_SEGMENT_LEN'],
        overlap = config['params']['NANOPOLISH']['NANOPOLISH_OVERLAP_LEN'],
        options = config['params']['NANOPOLISH']['OPTIONS'],
        liste_segments = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/NANOPOLISH/segments.txt",
        vcf = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/NANOPOLISH/polished.vcf",
        preset = f"{config['params']['MINIMAP2']['PRESET_OPTION']}" if {config['params']['MINIMAP2']['PRESET_OPTION']}!='' else 'map-pb'
    log:
        output = f"{output_dir}LOGS/CORRECTION/NANOPOLISH/{{fastq}}_{{assemblers}}_NANOPOLISH.o",
        error = f"{output_dir}LOGS/CORRECTION/NANOPOLISH/{{fastq}}_{{assemblers}}_NANOPOLISH.e",
    benchmark:
        f"{output_dir}LOGS/CORRECTION/NANOPOLISH/{{fastq}}_{{assemblers}}_NANOPOLISH-BENCHMARK.txt",
    priority: 30
    message:
        """
        Launching {rule}
        input:
            draft : {input.draft}
        output:
            fasta : {output.fasta}
        params:
            fastq : {params.fastq}
            fast5 : {params.fast5}
            directory : {params.directory}
            segment : {params.segment}
            overlap: {params.overlap}
            options: {params.options}
            liste segments : {params.liste_segments}
            vcf : {params.vcf}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['NANOPOLISH_SIMG']
    shell:
        """
        seqtk seq -A {params.fastq} > {output.temp_fasta}
        nanopolish index -d {params.fast5} {output.temp_fasta} 1>>{log.output} 2>>{log.error}
        minimap2 -ax {params.preset} -t {threads} {input.draft} {output.temp_fasta} | samtools sort -o {output.temp_bam} -T reads.tmp 1>>{log.output} 2>>{log.error}
        samtools index {output.temp_bam} 1>>{log.output} 2>>{log.error}
        python3 /nanopolish/scripts/nanopolish_makerange.py {input.draft} --segment-length {params.segment} --overlap-length {params.overlap} > {params.liste_segments} 2>> {log.error} #TODO ce chemin peut changer
        while read LINE; do echo "$LINE"; nanopolish variants --consensus -t {threads} -o {params.vcf}-"$LINE" -r {output.temp_fasta} -b {output.temp_bam} -g {input.draft} {params.options} -w "$LINE"; done < {params.liste_segments} 1>>{log.output} 2>>{log.error}
        nanopolish vcf2fasta --skip-checks -g {input.draft} {params.vcf}-* 1>{output.fasta} 2>>{log.error}
        """

rule run_medaka_train:
    """
    launching Medaka Train with fasta ref
    """
    threads: get_threads('run_medaka_train', 8)
    input:
        draft = draft_to_correction,
        fastq = get_fastq,
        ref = ref,
    output:
        fasta_cat_acc = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/MEDAKA/training/model.best.cat_acc.hdf5",
        fasta_val_cat_acc = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/MEDAKA/training/model.best.val_cat_acc.hdf5"
    params:
        index_fai = rules.index_fasta_to_correction.output.index_fai,
        index_fmmi = rules.index_fasta_to_correction.output.index_mmi,
        out_name = directory(f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/MEDAKA/"),
    log:
        output = f"{output_dir}LOGS/CORRECTION/MEDAKA/{{fastq}}_{{assemblers}}_MEDAKA_TRAIN.o",
        error = f"{output_dir}LOGS/CORRECTION/MEDAKA/{{fastq}}_{{assemblers}}_MEDAKA_TRAIN.e",
    benchmark:
        f"{output_dir}LOGS/CORRECTION/MEDAKA/{{fastq}}_{{assemblers}}_MEDAKA_TRAIN-BENCHMARK.txt",
    priority: 30
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
            draft : {input.draft}
            fastq : {input.fastq}
            ref: {input.ref}
        output:
            fasta_cat_acc: {output.fasta_cat_acc}
            fasta_val_cat_acc: {output.fasta_val_cat_acc}
        params:
            index_fai : {params.index_fai}
            index_fmmi : {params.index_fmmi}
            out_name: {params.out_name}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['MEDAKA_SIMG']
    shell:
        """
        mini_align -t {threads} -m -r {input.draft} -i {input.fastq} -p {params.out_name}calls2draft 1>{log.output} 2>{log.error}
        mini_align -t {threads} -m -r {input.draft} -i {input.ref} -p {params.out_name}truth2draft 1>>{log.output} 2>>{log.error}
        medaka features {params.out_name}calls2draft.bam {params.out_name}train_features.hdf --truth {params.out_name}truth2draft.bam --threads {threads} --batch_size 100 --chunk_len 1000 --chunk_ovlp 0 1>>{log.output} 2>>{log.error}
        medaka train {params.out_name}train_features.hdf --train_name {params.out_name}training --epochs 10 1>>{log.output} 2>>{log.error}
        """

rule run_medaka_consensus:
    """
    launching Medaka Consensus
    """
    threads: get_threads('run_medaka_consensus', 8)
    input:
        draft = draft_to_correction,
        fastq = get_fastq,
        model = f"{rules.run_medaka_train.output.fasta_val_cat_acc if config['params']['MEDAKA']['MEDAKA_TRAIN_WITH_REF'] else config['params']['MEDAKA']['MEDAKA_MODEL_PATH']}"
    output:
        fasta = f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/MEDAKA/consensus.fasta"
    params:
        dir = directory(f"{output_dir}{{fastq}}/{{assemblers}}/CORRECTION/MEDAKA"),
    log:
        output = f"{output_dir}LOGS/CORRECTION/MEDAKA/{{fastq}}_{{assemblers}}_MEDAKA_CONSENSUS.o",
        error = f"{output_dir}LOGS/CORRECTION/MEDAKA/{{fastq}}_{{assemblers}}_MEDAKA_CONSENSUS.e",
    benchmark:
        f"{output_dir}LOGS/CORRECTION/MEDAKA/{{fastq}}_{{assemblers}}_MEDAKA_CONSENSUS-BENCHMARK.txt"
    priority: 30
    message:
        """
        Launching {rule}.
        threads : {threads}
        input:
            draft : {input.draft}
            fastq : {input.fastq}
            model : {input.model}
        output:
            fasta : {output.fasta}
        params:
            dir : {params.dir}
        log:
            output : {log.output}
            error: {log.error}
        command :
        medaka_consensus -t {threads} -i {input.fastq} -d {input.draft} -o {params.dir} -m {input.model}
        """
    singularity:
        config['tools']['MEDAKA_SIMG']
    shell:
        """
        medaka_consensus -t {threads} -i {input.fastq} -d {input.draft} -o {params.dir} -m {input.model} 1>{log.output} 2>{log.error}
        """

def fasta_to_busco(wildcards):
    if wildcards.busco_step == 'STEP_ASSEMBLY':
        return f"{output_dir}{{fastq}}/{{assemblers}}/ASSEMBLER/assembly{add_circular_name}.fasta"
    elif wildcards.busco_step == 'STEP_POLISHING_RACON':
        return f"{output_dir}{{fastq}}/{{assemblers}}/POLISHING/RACON/racon_{nb}/assembly.racon{nb}.fasta"
    elif wildcards.busco_step == 'STEP_CORRECTION_NANOPOLISH':
        return rules.run_nanopolish.output.fasta
    elif wildcards.busco_step == 'STEP_CORRECTION_MEDAKA':
        return rules.run_medaka_consensus.output.fasta
    else:
        raise ValueError("problem with fasta to busco rule.")


################################ QUALITY CHECK ####################################

rule preparing_fasta_to_quality:
    """
    preparing fasta to quality
    """
    threads: get_threads('preparing_fasta_to_quality', 2)
    input:
        fasta = fasta_to_busco
    output:
        renamed = f"{output_dir}{{fastq}}/QUAST/data/{{assemblers}}-{{busco_step}}-{add_circular_name}assembly.fasta",
    params:
        #index_fai = f"{output_dir}{{fastq}}/QUAST/data/{{assemblers}}-{{busco_step}}-{add_circular_name}assembly.fasta.fai",
        index_mmi = f"{output_dir}{{fastq}}/QUAST/data/{{assemblers}}-{{busco_step}}-{add_circular_name}assembly.fasta.mmi"
    log:
        output = f"{output_dir}LOGS/QUALITY/QUAST/{{fastq}}-{{assemblers}}-{{busco_step}}-preparing_QUAST.o",
        error = f"{output_dir}LOGS/QUALITY/QUAST/{{fastq}}-{{assemblers}}-{{busco_step}}-preparing_QUAST.e",
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fasta = {input.fasta},
        output:
            renamed : {output.renamed}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['MINIMAP2_SIMG']
    shell:
        """
        ln -s {input.fasta} {output.renamed} 1>{log.output} 2>{log.error}
        samtools faidx {output.renamed} 1>>{log.output} 2>>{log.error}
        minimap2 -d {params.index_mmi} {output.renamed} 1>>{log.output} 2>>{log.error}
        """

rule run_quast:
    """
    preparing fasta to quast and launch quast
    """
    threads: get_threads('run_quast', 4)
    input:
        liste = expand(rules.preparing_fasta_to_quality.output.renamed, fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=BUSCO_STEPS),
    output:
        report = f"{output_dir}{{fastq}}/QUAST/quast_results/report.html",
        report_path = f"{output_dir}REPORT/{{fastq}}/report_quast.html"
    params:
        genome_size = f"{config['params']['QUAST']['GENOME_SIZE_PB']}",
        #indir = f"{output_dir}{{fastq}}/QUAST/data/",
        directory = f"{output_dir}{{fastq}}/QUAST/quast_results/",
        reference = f"{config['DATA']['REF']}" if f"{config['DATA']['REF']}" == '' else f" -r {config['DATA']['REF']}",
        gff = f"{config['params']['QUAST']['GFF']}" if f"{config['params']['QUAST']['GFF']}"=='' else f" --g {config['params']['QUAST']['GFF']}",
        options = f"{config['params']['QUAST']['OPTIONS']}",
    log:
        output = f"{output_dir}LOGS/QUALITY/QUAST/QUAST_{{fastq}}.o",
        error = f"{output_dir}LOGS/QUALITY/QUAST/QUAST_{{fastq}}.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/QUAST/QUAST_{{fastq}}-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            liste = {input.liste},
        output:
            report: {output.report}
        params:
            genome_size  : {params.genome_size}
            directory : {params.directory}
            reference = {params.reference},
            gff : {params.gff}
            options : {params.options}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['QUAST_SIMG']
    shell:
        """
        quast.py -t {threads} -o {params.directory} {params.gff} {params.reference} {params.options} {input.liste} 1>{log.output} 2>{log.error};
        cp {output.report} {output.report_path}
        """

rule run_busco:
    """
    BUSCO v4 assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs v10
    """
    threads: get_threads('run_busco', 4)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed
    output:
        summary = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/BUSCO_RESULTS/short_summary_BUSCO_RESULTS.txt",
    params:
        out_path = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/",
        busco_name = f"BUSCO_RESULTS",
        database = config['params']['BUSCO']['DATABASE'],
        model = config['params']['BUSCO']['MODEL'],
        sp = f" -sp {config['params']['BUSCO']['SP']}" if f"{config['params']['BUSCO']['SP']}"!="" else f"{config['params']['BUSCO']['SP']}",
        #path_tool = config['tools']['BUSCO_TOOL']
    log:
        output = f"{output_dir}LOGS/QUALITY/BUSCO/{{fastq}}_{{assemblers}}_{{busco_step}}_BUSCO.o",
        error = f"{output_dir}LOGS/QUALITY/BUSCO/{{fastq}}_{{assemblers}}_{{busco_step}}_BUSCO.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/BUSCO/{{fastq}}_{{assemblers}}_{{busco_step}}_BUSCO-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule}
        threads : {threads}
        input:
            fasta : {input.fasta}
        output:
            summary : {output.summary}
        params:
            out_path : {params.out_path}
            busco_name : {params.busco_name}
            database : {params.database}
            model : {params.model}
            sp : {params.sp}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['BUSCO_SIMG']
    shell:
        """
        busco -f -i {input.fasta} -o {params.busco_name} --out_path {params.out_path} -l {params.database} --offline -m {params.model} {params.sp} -c {threads} 1>{log.output} 2>{log.error}
        mv {params.out_path}BUSCO_RESULTS/short_summary*BUSCO_RESULTS.txt {params.out_path}BUSCO_RESULTS/short_summary_BUSCO_RESULTS.txt 1>>{log.output} 2>>{log.error}
        """

rule run_diamond:
    """
    running diamond to blobtools
    """
    threads: get_threads('run_diamond', 4)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed,
    output:
        csv = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/DIAMOND/diamond.csv",
    params:
        db = f"{config['params']['DIAMOND']['DATABASE']}",
        dir = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/DIAMOND/",
        mode = "blastx",
        format = '6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        #path_tool = config['tools']['DIAMOND_TOOL']
    log:
        output = f"{output_dir}LOGS/QUALITY/DIAMOND/{{fastq}}_{{assemblers}}_{{busco_step}}_DIAMOND.o",
        error = f"{output_dir}LOGS/QUALITY/DIAMOND/{{fastq}}_{{assemblers}}_{{busco_step}}_DIAMOND.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/DIAMOND/{{fastq}}_{{assemblers}}_{{busco_step}}_DIAMOND-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fasta : {input.fasta}
        output:
            csv : {output.csv}
        params:
            database: {params.db}
            mode: {params.mode}
            format : {params.format}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['DIAMOND_SIMG']
    shell:
        """
        mkdir -p {params.dir}
        diamond {params.mode} --query {input.fasta} --db {params.db} --outfmt {params.format} --threads {threads} --sensitive --max-target-seqs 1 --evalue 1e-25 --threads {threads} --out {output.csv} 1>{log.output} 2>{log.error}
        """

rule run_minimap2:
    """
    running minimap2 in mode mapping ONT
    """
    threads: get_threads('run_minimap2', 8)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed,
        fastq = get_fastq,
    output:
        bam = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/MINIMAP2/minimap2mapping.bam",
    params:
        preset = f"{config['params']['MINIMAP2']['PRESET_OPTION']}" if {config['params']['MINIMAP2']['PRESET_OPTION']}!='' else 'map-pb',
    log:
        output = f"{output_dir}LOGS/QUALITY/MINIMAP2/{{fastq}}_{{assemblers}}_{{busco_step}}_MINIMAP2.o",
        error = f"{output_dir}LOGS/QUALITY/MINIMAP2/{{fastq}}_{{assemblers}}_{{busco_step}}_MINIMAP2.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/MINIMAP2/{{fastq}}_{{assemblers}}_{{busco_step}}_MINIMAP2-BENCHMARK.txt"
    #params:
    #    index_fai = rules.preparing_fasta_to_quality.output.index_fai,
    #    index_mmi = rules.preparing_fasta_to_quality.output.index_mmi,
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fasta : {input.fasta}
            fastq : {input.fastq}
        output:
            bam : {output.bam}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['NANOPOLISH_SIMG']
    shell: #TODO: gerer la sortie standard de minimap2 dans les logs
        """
        minimap2 -x {params.preset} -t {threads} {input.fasta} {input.fastq} -a | samtools sort -o {output.bam} 1>{log.output} 2>{log.error}
        samtools index {output.bam}  1>>{log.output} 2>>{log.error}
        """

rule run_blobtools:
    """
    blobtools v1
    """
    threads: get_threads('run_blobtools', 8)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed,
        sorted_bam = rules.run_minimap2.output.bam,
        diamond = rules.run_diamond.output.csv,
    output:
        table = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/BLOBTOOLS/output.quality.blobDB.table.txt",
        read_cov = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/BLOBTOOLS/read_cov.png",
        blob = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/BLOBTOOLS/blob.png"
    params:
        dir = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/BLOBTOOLS/",
        #TODO mettre en config
        names = '/usr/local/blobtools-blobtools_v1.1.1/data/names.dmp',
        nodes = '/usr/local/blobtools-blobtools_v1.1.1/data/nodes.dmp',
    log:
        output = f"{output_dir}LOGS/QUALITY/BLOBTOOLS/{{fastq}}_{{assemblers}}_{{busco_step}}_BLOBTOOLS.o",
        error = f"{output_dir}LOGS/QUALITY/BLOBTOOLS/{{fastq}}_{{assemblers}}_{{busco_step}}_BLOBTOOLS.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/BLOBTOOLS/{{fastq}}_{{assemblers}}_{{busco_step}}_BLOBTOOLS-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fasta : {input.fasta}
            sorted_bam : {input.sorted_bam}
            diamond : {input.diamond}
        output:
            table: {output.table}
        params:
            names : {params.names}
            nodes : {params.nodes}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['BLOBTOOLS_SIMG']
    shell:
        """
        cd {params.dir} 1>{log.output} 2>{log.error}
        blobtools create -i {input.fasta} -b {input.sorted_bam} -t {input.diamond} -o quality --names {params.names} --nodes {params.nodes}  1>>{log.output} 2>>{log.error}
        blobtools view -i quality.blobDB.json --cov -o output 1>>{log.output} 2>>{log.error}
        blobtools plot -i quality.blobDB.json 1>>{log.output} 2>>{log.error}
        mv quality.blobDB.*.blobplot.read_cov.bam0.png {output.read_cov} 1>>{log.output} 2>>{log.error}
        mv quality.blobDB.*.blobplot.bam0.png {output.blob} 1>>{log.output} 2>>{log.error}
        """

rule run_weesam:
    """
    weesam runs only in the last assemblies
    """
    threads: get_threads('run_weesam', 4)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed,
        fastq = get_fastq,
        bam = rules.run_minimap2.output.bam,
    output:
        txt = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/WEESAM/minimap2mapping.txt",
        weesam_outdir = directory(f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/WEESAM/minimap2mapping_html_results")
    params:
        out_dir = lambda w, output: os.path.dirname(output.txt),
    log:
        output = f"{output_dir}LOGS/QUALITY/WEESAM/{{fastq}}_{{assemblers}}_{{busco_step}}_WEESAM.o",
        error = f"{output_dir}LOGS/QUALITY/WEESAM/{{fastq}}_{{assemblers}}_{{busco_step}}_WEESAM.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/WEESAM/{{fastq}}_{{assemblers}}_{{busco_step}}_WEESAM-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fastq : {input.fastq}
            fasta : {input.fasta}
            bam : {input.bam}
        output:
            txt: {output.txt}
        params:
            out_dir : {params.out_dir}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['WEESAM_SIMG']
    shell:
        """
        weeSAM --overwrite --bam {input.bam} --html "{params.out_dir}"/minimap2mapping.html \
            --out "{params.out_dir}"/minimap2mapping.txt 1>>{log.output} 2>>{log.error}
        """

rule run_mummer:
    """
    This rule run nucmer for assemblytics
    """
    threads: get_threads('run_mummer', 4)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed,
        ref = ref
    output:
        delta = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/MUMMER/mummer.delta.gz",
    params:
        minmatch = config['params']['MUMMER']['MINMATCH'], #option -l  20
        mincluster = config['params']['MUMMER']['MINCLUSTER'], #option -c 65
        prefix = "mummer",
        dir = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/MUMMER/",
 #       path_tool = config['tools']['MUMMER_TOOL']
    log:
        output = f"{output_dir}LOGS/QUALITY/MUMMER/{{fastq}}_{{assemblers}}_{{busco_step}}_MUMMER.o",
        error = f"{output_dir}LOGS/QUALITY/MUMMER/{{fastq}}_{{assemblers}}_{{busco_step}}_MUMMER.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/MUMMER/{{fastq}}_{{assemblers}}_{{busco_step}}_MUMMER-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fasta : {input.fasta}
            ref : {input.ref}
        output:
            delta: {output.delta}
        params:
            minmatch = {params.minmatch}
            mincluster = {params.mincluster}
            prefix: {params.prefix}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['MUMMER_SIMG']
    shell:
        """
        cd {params.dir} 1>{log.output} 2>{log.error}
        nucmer --maxmatch  -l {params.minmatch} -c {params.mincluster} {input.ref} {input.fasta} -p {params.prefix} 1>>{log.output} 2>>{log.error}
        gzip mummer.delta 1>>{log.output} 2>>{log.error}
        """

rule run_assemblytics:
    """
    Assemblytics analyze your assembly by comparing it to a reference genome https://github.com/MariaNattestad/assemblytics
    """
    threads: get_threads('run_assemblytics', 4)
    input:
        delta = rules.run_mummer.output.delta,
    output:
        summary = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/ASSEMBLYTICS/OUT.Assemblytics_structural_variants.summary",
        png_dotplot = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/ASSEMBLYTICS/OUT.Assemblytics.Dotplot_filtered.png",
        png_Nchart = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/ASSEMBLYTICS/OUT.Assemblytics.Nchart.png",
        png_log_all_sizes = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/ASSEMBLYTICS/OUT.Assemblytics.size_distributions.all_variants.log_all_sizes.png",
    params:
        dir = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/ASSEMBLYTICS/",
        unique_anchor_length = config['params']['ASSEMBLYTICS']['UNIQUE_ANCHOR_LEN'],
        min_variant_size = config['params']['ASSEMBLYTICS']['MIN_VARIANT_SIZE'],
        max_variant_size = config['params']['ASSEMBLYTICS']['MAX_VARIANT_SIZE'],
        prefix = "OUT",
    log:
        output = f"{output_dir}LOGS/QUALITY/ASSEMBLYTICS/{{fastq}}_{{assemblers}}_{{busco_step}}_ASSEMBLYTICS.o",
        error = f"{output_dir}LOGS/QUALITY/ASSEMBLYTICS/{{fastq}}_{{assemblers}}_{{busco_step}}_ASSEMBLYTICS.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/ASSEMBLYTICS/{{fastq}}_{{assemblers}}_{{busco_step}}_ASSEMBLYTICS-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            delta : {input.delta}
        output:
            summary: {output.summary}
        params:
            unique_anchor_length = {params.unique_anchor_length}
            min_variant_size = {params.min_variant_size}
            max_variant_size = {params.max_variant_size}
            prefix: {params.prefix}
        log:
            output : {log.output}
            error: {log.error}
        """
    singularity:
        config['tools']['ASSEMBLYTICS_SIMG']
    shell:
        """
        cd {params.dir} 1>{log.output} 2>{log.error}
        Assemblytics {input.delta} {params.prefix} {params.unique_anchor_length} {params.min_variant_size} {params.max_variant_size} 1>>{log.output} 2>>{log.error}
        """

rule run_flagstat:
    """
    calculate stats from mapping: use to quality report
    """
    threads: get_threads('run_flagstat', 4)
    input:
        bam = rules.run_minimap2.output.bam,
    output:
        txt = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/FLAGSTAT/flagstat.txt",
    log:
        error = f"{output_dir}LOGS/QUALITY/FLAGSTAT/{{fastq}}_{{assemblers}}_{{busco_step}}_FLAGSTAT.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/FLAGSTAT/{{fastq}}_{{assemblers}}_{{busco_step}}_FLAGSTAT-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            bam : {input.bam}
        output:
            txt: {output.txt}
        log:
            error: {log.error}
        """
    singularity:
        config['tools']['MEDAKA_SIMG']
    shell:
        """
        samtools flagstat {input.bam} --threads {threads} > {output.txt} 2>{log.error}
        """

rule combined_fastq:
    """
    zcat des fastq.gz. fastq sequences must to be decompressed to KAT.
    """
    threads: get_threads('combined_fastq', 2)
    input:
        illumina_rep = illumina
    params:
        command = "zcat " if ext_illumina == '.gz' else  "cat "
    output:
        combined_data = f"{output_dir}combined_data.fastq",
    log:
        error = f"{output_dir}LOGS/QUALITY/KAT/combinedfastq2KAT.e",
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            illumina_rep :{input.illumina_rep}
        output:
            combined: {output.combined_data}
        log:
            error: {log.error}
        command:
            {params.command} {input.illumina_rep}*{ext_illumina} > {output.combined_data}
        """
    shell:
        """
            {params.command} {input.illumina_rep}*{ext_illumina} > {output.combined_data}
        """

rule run_KAT:
    """
    KAT is useful tool for high accuracy sequence data.
    The spectra-cn (copy number spectra) graph shows a decomposition of k-mers in the assembly vs k-mers in the reads.
    The black portion are k-mers not present in the assembly, the red portion is found once in the assembly, and so on.
    This shows the completeness of an assembly, i.e. are all the reads assembled into contigs representative of the sequence data.
    https://kat.readthedocs.io/en/latest/using.html
    """
    threads: get_threads('run_KAT', 4)
    input:
        fasta = rules.preparing_fasta_to_quality.output.renamed,
        combined_data = rules.combined_fastq.output.combined_data,
    output:
        gcp = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.gcp", #matrix
        png_hist = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.hist.png",
        png_spectra_hist = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat-spectra-hist.png",
        png_spectra_cn = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat-spectra-cn.png",
        png_density = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat-density.png",
        png_gcp_mx = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.gcp.mx.png",
        png_comp_main_mx_spectra_cn = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.comp-main.mx.spectra-cn.png",
    params:
        hist = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.hist",
        cmp = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.comp", #matrix
        sect = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.sect",
        dir = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/",
        gcp_tmp = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.gcp.mx", #matrix
        cmp_tmp  = f"{output_dir}{{fastq}}/{{assemblers}}/QUALITY/{{busco_step}}/KAT/kat.comp-main.mx", #matrixkat.comp-main.mx
        #path_tool = config['tools']['KAT_TOOL'],
    log:
        output = f"{output_dir}LOGS/QUALITY/KAT/{{fastq}}_{{assemblers}}_{{busco_step}}_KAT.o",
        error = f"{output_dir}LOGS/QUALITY/KAT/{{fastq}}_{{assemblers}}_{{busco_step}}_KAT.e",
    benchmark:
        f"{output_dir}LOGS/QUALITY/KAT/{{fastq}}_{{assemblers}}_{{busco_step}}_KAT-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            fasta  : {input.fasta}
            combined_data : {input.combined_data}
        params:
            dir : {params.dir}
            cmp : {params.cmp}
            hist : {params.hist}
            hist : {params.sect}
        output:
            cmp : {output.gcp}
        log:
            output : {log.output}
            error : {log.error}
        """
    singularity:
        config['tools']['KAT_SIMG']
    shell:
        """
        cd {params.dir} 1>{log.output} 2>{log.error}
        kat hist -t {threads} -o {params.hist} {input.combined_data} 1>{log.output} 2>{log.error}
        kat gcp -t {threads} -o {output.gcp} {input.combined_data} 1>>{log.output} 2>>{log.error}
        kat comp -t {threads} -o {params.cmp} {input.combined_data} {input.fasta} 1>>{log.output} 2>>{log.error}
        kat plot spectra-hist {params.hist} 1>>{log.output} 2>>{log.error}
        kat plot density {params.gcp_tmp} 1>>{log.output} 2>>{log.error}
        kat plot density {params.cmp_tmp} 1>>{log.output} 2>>{log.error}
        kat plot spectra-cn {params.cmp_tmp} 1>>{log.output} 2>>{log.error}
        mv {params.gcp_tmp} {output.gcp} 1>>{log.output} 2>>{log.error}
        """

######### Standardizing starting coordinate of bacterial genomes ############

rule run_fixstart:
    """
    Standardizing starting coordinate of bacterial genome assemblies with fixstart module of Circlator.
    """
    threads: get_threads('run_fixstart', 1)
    input:
        assembly_file = rules.preparing_fasta_to_quality.output.renamed
    output:
        fix_start_fasta = f"{output_dir}{{fastq}}/{{assemblers}}/MSA/FIXSTART-{{busco_step}}/startfixed_asm.fasta",
    params:
        fix_start_log = f"{output_dir}{{fastq}}/{{assemblers}}/MSA/FIXSTART-{{busco_step}}/startfixed_asm.log",
        out_dir = lambda w, output: os.path.dirname(output.fix_start_fasta)
    log:
        output = f"{output_dir}LOGS/FIXSTART/{{fastq}}_{{assemblers}}_{{busco_step}}-FIXSTART.o",
        error =  f"{output_dir}LOGS/FIXSTART/{{fastq}}_{{assemblers}}_{{busco_step}}-FIXSTART.e",
    benchmark:
        f"{output_dir}LOGS/FIXSTART/{{fastq}}_{{assemblers}}_{{busco_step}}-FIXSTART-BENCHMARK.txt"
    priority: 20
    message:
        """
        Launching {rule} ...
        threads : {threads}
        input:
            assembly_file  : {input.assembly_file}
        params:
            fix_start_log : {params.fix_start_log}
            out_dir : {params.out_dir}
        output:
            fix_start_fast : {output.fix_start_fasta}
        """
    singularity:
        config['tools']['MINICONDA_SIMG']
    conda:
         "envs/run_circlator_cenv.yml"
    shell:
        """
        exec > >(tee "{log.output}") 2>&1
        
        # inspect fasta titles for circular molecules and write their names in a file
        list_lin_seqs()
        {{
          grep -o -E "^>.*$" $1 | grep -v -E "circular" | tr -d ">" > "$2"
        }}
        linSeqNamesFile="{params.out_dir}/linSeqNames.txt"
        set +e
        list_lin_seqs "{input.assembly_file}" "$linSeqNamesFile"
        set -e
        echo "##  $(date): processing {input.assembly_file} "
        echo "##  Using circlator fixstart to rotate circular sequences so that they start at a dnaA gene (if found)"
        echo "##  Ignoring the following linear sequences: $(cat $linSeqNamesFile)."

        circlator fixstart --ignore="$linSeqNamesFile" "{input.assembly_file}" "{params.out_dir}"/startfixed_asm

        exit 0
        """

rule run_mauve:
    threads:
        get_threads('run_mauve', 8)
    input:
        #liste = expand(f"{rules.run_fixstart.output.fix_start_fasta}", fastq = FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step=input_last()) if config['MSA']['FIXSTART'] else expand(f"{rules.preparing_fasta_to_quality.output.renamed}", fastq = FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step=input_last())
        liste = expand(f"{rules.run_fixstart.output.fix_start_fasta}", fastq = FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step=BUSCO_STEPS) if config['MSA']['FIXSTART'] else expand(f"{rules.preparing_fasta_to_quality.output.renamed}", fastq = FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step=BUSCO_STEPS)

    params:
        out_dir = f"{output_dir}{{fastq}}/MAUVE_ALIGN/",
        circular = config['DATA']['CIRCULAR'],
        dir_quast = f"{output_dir}{{fastq}}/QUAST/data/"
    output:
        xmfa = f"{output_dir}{{fastq}}/MAUVE_ALIGN/candidate_assemblies.xmfa",
    priority: 20
    message:
        """
        Launching {rule}
        Generating Mauve multiple alignment of sequences in all final assembly fasta files
        threads : {threads}
        input:
            liste: {input.liste}
        output:
            xmfa: {output.xmfa}
        """
    log:
        output = f"{output_dir}LOGS/MAUVE/{{fastq}}-MAUVE.o",
        error = f"{output_dir}LOGS/MAUVE/{{fastq}}-MAUVE.e",
    singularity:
        config['tools']['MINICONDA_SIMG']
    conda:
        'envs/mauve_cenv.yml'
    script:
         "AdditionalScripts/run_mauve.py"

    #shell:
    #    """
    #    progressiveMauve --output={output.xmfa} {input.liste}
    #    """

################################ RAPPORT ####################################

rule rule_graph:
    """
    run dag on {rule}
    """
    threads: get_threads('rule_graph', 1)
    input:
        conf = str(path_config),
    params:
        tmp = f"{output_dir}dag.tmp"
    output:
        dag = f"{output_dir}dag.png"
    log:
        output = f"{output_dir}LOGS/REPORT/dag.o",
        error = f"{output_dir}LOGS/REPORT/dag.e",
    priority: 10
    message:
        """
        making dag ...
        snakemake -n --rulegraph --configfile {input.conf} > {params.tmp} 2>{log.error}
        dot -Tpng {params.tmp} >{output.dag} 2>>{log.error}
        """
    singularity:
        config['tools']['MINICONDA_SIMG']
    conda:
        "envs/environment.yaml"
    shell:
        """
        snakemake -n --rulegraph --configfile {input.conf} >{params.tmp} 2>{log.error}
        dot -Tpng {params.tmp} >{output.dag} 2>>{log.error}
        rm {params.tmp} 2>>{log.output} 2>>{log.error}
        """

rule run_stats:
    """
    run stat
    """
    threads: get_threads('run_stats', 1)
    input:
        fastq = get_fastq,
        summary = expand(rules.run_busco.output.summary, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = BUSCO_STEPS)
    output:
        stat = f"{output_dir}REPORT/{{fastq}}/Stats_busco.csv",
    params:
        liste_busco = expand(BUSCO_STEPS),
        liste_assemblers = expand(ASSEMBLY_TOOLS),
        out_dir = directory(f"{output_dir}{{fastq}}"),
    priority: 10
    message:
        """
        make stats
        """
    script:
        "reports/getStats.py"

rule run_benchmark_time:
    """
    run benchmark
    """
    threads: get_threads('run_benchmark_time', 1)
    input:
        fastq = get_fastq,
        summary = expand(rules.run_busco.output.summary, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = BUSCO_STEPS)
    output:
        stat = f"{output_dir}REPORT/{{fastq}}/Stats_benchmark.csv",
    params:
        liste_busco = expand(BUSCO_STEPS),
        liste_assemblers = expand(ASSEMBLY_TOOLS),
        out_dir = directory(f"{output_dir}"),
        racon_round = config['params']['RACON']['RACON_ROUNDS'],
        liste_correcteurs = expand(CORRECTION_TOOLS),
        medaka_training = config['params']['MEDAKA']['MEDAKA_TRAIN_WITH_REF']
    priority: 10
    message:
        """
        make benchmark
        """
    script:
        "reports/getBenchmark.py"

#TODO : quast could be optional
rule run_report:
    """
    print report
    """
    threads: get_threads('run_report', 1)
    input:
        stat =  expand(rules.run_stats.output.stat,fastq=FASTQ),
        bench = expand(rules.run_benchmark_time.output.stat,fastq=FASTQ),
        quast =  expand(rules.run_quast.output.report_path,fastq=FASTQ),
        dag = rules.rule_graph.output.dag,
        check = rules.run_dico_final.output.check
    params:
        #dag = rules.rule_graph.output.dag,
        conf = path_config,
        outdir = directory(f"{output_dir}"),
        out_dir = directory(f"{output_dir}REPORT"),
        liste_assembler = ASSEMBLY_TOOLS,
        #list_final = output_final,
        weesam_outdir = expand(rules.run_weesam.output.weesam_outdir , fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['WEESAM'] else '',
        blob = expand(rules.run_blobtools.output.blob, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['BLOBTOOLS'] else '',
        blob_read_cov = expand(rules.run_blobtools.output.read_cov, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['BLOBTOOLS'] else '',
        weesam_txt = expand(rules.run_weesam.output.txt, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['WEESAM'] else '',
        kat_png_hist = expand(rules.run_KAT.output.png_hist, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['KAT'] else '',
        kat_png_spectra_hist = expand(rules.run_KAT.output.png_spectra_hist, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['KAT'] else '',
        kat_png_spectra_cn = expand(rules.run_KAT.output.png_spectra_cn, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['KAT'] else '',
        kat_png_density = expand(rules.run_KAT.output.png_density, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['KAT'] else '',
        kat_png_gcp_mx = expand(rules.run_KAT.output.png_gcp_mx, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['KAT'] else '',
        kat_png_comp_main_mx_spectra_cn = expand(rules.run_KAT.output.png_comp_main_mx_spectra_cn, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['KAT'] else '',
        assemblytics_png_dotplot = expand(rules.run_assemblytics.output.png_dotplot, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['ASSEMBLYTICS'] else '',
        assemblytics_png_Nchart = expand(rules.run_assemblytics.output.png_Nchart, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['ASSEMBLYTICS'] else '',
        assemblytics_png_log_all_sizes = expand(rules.run_assemblytics.output.png_log_all_sizes, fastq=FASTQ, assemblers = ASSEMBLY_TOOLS, busco_step = input_last()) if config['QUALITY']['ASSEMBLYTICS'] else '',
        mauve = expand(rules.run_mauve.output.xmfa, fastq=FASTQ) if config['MSA']['MAUVE'] else '',
        #fixstart = expand(rules.run_fixstart.output.fix_start_fasta, fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last()) if (config['MSA']['FIXSTART'] and not config['MSA']['MAUVE']) else '',
        ##fixstart = expand(rules.run_fixstart.output.fix_start_fasta, fastq=FASTQ, assemblers=ASSEMBLY_TOOLS, busco_step=input_last()) if (config['DATA']['CIRCULAR'] and config['MSA']['FIXSTART'] and not config['MSA']['MAUVE']) else '',
    output:
        report = f"{output_dir}REPORT/Report.html",
    #priority: 0
    message:
        """
        print final report...
        """
    singularity:
        config['tools']['R_SIMG']
    script:
        "reports/Report.Rmd"
