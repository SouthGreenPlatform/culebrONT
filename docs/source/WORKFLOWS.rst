How to create a workflow ?
==========================

CulebrONT allow you to build a workflow using a simple configuration ``config.yaml`` file. In this file :

* First, provide data paths
* Second, activate tools from assembly to correction.
* Third, activate tools from quality checking of assemblies.
* And last, manage parameters tools.

1. Providing data
------------------

First, indicate the data path in the configuration ``config.yaml`` file:

.. code-block:: yaml

   DATA:
       FASTQ: '/path/to/fastq/directory/'
       REF: '/path/to/referencefile.fasta'
       GENOME_SIZE: '1m'
       FAST5: '/path/to/fast5/directory/'
       ILLUMINA: '/path/to/illumina/directory/'
       OUTPUT: '/path/to/output/directory/'

Find here a summary table with description of each data need to lauch CulebrONT :

.. csv-table::
    :header: "Input", "Description"
    :widths: auto

    "FASTQ", "Every FASTQ file should contain the whole set of reads to be assembled. Each fastq file will be assembled independently"
    "REF","Only one REFERENCE genome file will be used by CulebrONT. This REFERENCE will be used for quality steps (ASSEMBLYTICS, QUAST and MAUVE)"
    "GENOME_SIZE", "Estimated genome size of the assembly can be done on mega (Mb), giga(Gb) or kilobases (Kb). This size is used on some assemblers (CANU) and also on QUAST quality step"
    "FAST5","Nanopolish needs FAST5 files to training steps. Please give the path of FAST5 folder in the *FAST5* DATA parameter. Inside this directory, a subdirectory with the exact same name as the corresponding FASTQ (before the *.fastq.gz*\ ) is requested. For instance, if in the *FASTQ* directory we have *run1.fastq.gz* and *run2.fastq.gz*\ , CulebrONT is expecting the *run1/* and *run2/* subdirectories in the FAST5 main directory"
    "ILLUMINA","Indicate the path to the directory with *Illumina* sequence data (in fastq or fastq.gz format) to perform KAT quality. Use preferentially paired-end data. All fastq files need to be homogeneous on their extension name"
    "OUTPUT","output *path* directory"

.. warning::

    For FASTQ, naming convention accepted by CulebrONT is *NAME.fastq.gz* or *NAME.fq.gz* or *NAME.fastq* or *NAME.fq*

    All fastq files have to be homogeneous on their extension and can be compressed or not.

    Reference fasta file need a fasta or fa extension uncompressed.

2. Chose assemblers, polisher and correctors
--------------------------------------------

Activate/deactivate assemblers, polishers and correctors as you wish.
Feel free to activate only assembly, assembly+polishing or assembly+polishing+correction.

.. note::
    If you are working on prokaryote, is recommendated to activate CIRCULAR steps

Example:

.. code-block:: yaml

   ASSEMBLY:
       CANU : False
       FLYE : True
       MINIASM : False
       SHASTA : False
       SMARDENOVO : True
       RAVEN: True
   CIRCULAR: False
   POLISHING:
       RACON: True
   CORRECTION:
       NANOPOLISH : True
       MEDAKA : False
   FIXSTART: False


3. Chose quality tools
----------------------
With CulebrONT you can use several quality tools to check assemblies.

* If BUSCO or QUAST are used, every fasta generated on the pipeline will be used with them.
* If BLOBTOOLS, ASSEMBLYTICS, WEESAM and KAT are activated only the last draft generated on the pipeline will be used.
* KAT quality tool can be activate but Illumina reads are mandatory in this case. These reads can be compressed or not.

.. code-block:: yaml

    QUALITY:
       BUSCO: True
       QUAST: True
       WEESAM: True
       BLOBTOOLS: True
       ASSEMBLYTICS: True
       KAT: True


Alignment of various assemblies **for small genomes (<10-20Mbp)** is also possible by using Mauve.

* If you want to improve alignment with MAUVE on circular molecules is recommended to activate *Fixstart* step.
* Only activate MAUVE if you have more than one sample and more than one quality step.

.. code-block:: yaml

   MSA:
       MAUVE: True


4. Parameters for some specific tools
--------------------------------------

You can manage tools parameters on the params section on ``config.yaml`` file.

Specifically to ``Racon``:

* Racon can be launch recursively from 1 to 9 rounds. 2 or 3 are recommended.

Specifically to ``Medaka``:

* If 'MEDAKA_TRAIN_WITH_REF' is activated, Medaka launchs training using the reference found in 'DATA/REF' param. Medaka does not take into account other medaka model parameters.

* If 'MEDAKA_TRAIN_WITH_REF' is deactivated, Medaka does not launch training but uses instead the model provided in 'MEDAKA_MODEL_PATH'. Give to CulebrONT path of medaka model in order to correct assemblies. This parameter could not be empty.

.. important::
    Medaka models can be downloaded from the medaka repository. You need to install ``git lfs`` (see documentation here https://git-lfs.github.com/) to download largest files before ``git clone https://github.com/nanoporetech/medaka.git\``.


Here you find standard parameters used on CulebrONT. Feel free to adapt it to your requires.

.. code-block:: yaml

    params:
        MINIMAP2:
            PRESET_OPTION: 'map-ont'
        FLYE:
            OPTIONS: ''
        CANU:
            MAX_MEMORY: '15G'
            OPTIONS: '-fast'
        SMARTDENOVO:
            KMER_SIZE: 16
            OPTIONS: '-J 5000'
        SHASTA:
            MEM_MODE: 'filesystem'
            MEM_BACKING: 'disk'
        CIRCLATOR:
            OPTIONS: ''
        RACON:
            RACON_ROUNDS: 2
        CORRECTION_MAKERANGE:
            SEGMENT_LEN: '50000'
            OVERLAP_LEN: '200'
        NANOPOLISH:
            OPTIONS: ''
        MEDAKA:
            MEDAKA_TRAIN_WITH_REF: True
            MEDAKA_MODEL_PATH: 'Data-Xoo-sub/medaka-models/r941_min_high_g303_model.hdf5'
            MEDAKA_FEATURES_OPTIONS: '--batch_size 100 --chunk_len 10000 --chunk_ovlp 1000'
            MEDAKA_TRAIN_OPTIONS: '--batch_size 100 --epochs 5000 '
            MEDAKA_CONSENSUS_OPTIONS: '-batch 50'
        BUSCO:
            DATABASE : 'Data-Xoo-sub/bacteria_odb10'
            MODEL : 'genome'
            SP : ''
        QUAST:
            GFF: ''
            OPTIONS : ''
        DIAMOND:
            DATABASE: 'Data-Xoo-sub/testBacteria.dmnd'
        MUMMER:
            MINMATCH : 100
            MINCLUSTER: 500
        ASSEMBLYTICS:
            UNIQUE_ANCHOR_LEN: 10000
            MIN_VARIANT_SIZE: 50
            MAX_VARIANT_SIZE: 10000

.. warning::
    Please check documentation of each tool and make sure that the settings are correct!

.. ############################################################

How to run the workflow ?
=========================

Command line
------------

Before lauch culebrONT please be sure you have already modified the ``config.yaml`` file as was explained on :ref:`WORKFLOWS:1. Providing data`

This is the recommended Snakemake command line to run CulebrONT:

.. code-block:: bash

    snakemake --nolock --use-conda --use-singularity --singularity-args '--bind $HOME' --cores -p -s Snakefile --latency-wait 6000000 --keep-going --restart-times 0 --rerun-incomplete --configfile config.yaml --conda-prefix $PWD/build_conda_envs

``config.yaml`` file is give to Snakemake by the argument ``--configfile``

To launch CulebrONT, you should use the parameters ``--use-singularity`` and ``--use-conda``. Please don't forget to export conda on your $PATH.

Snakemake compiles in each output directory conda environment. To avoid this, please use ``--conda-prefix /path/to/build_conda_env`` on snakemake command line.

Bind mount disks to singularity environment by using ``--singularity-args '--bind $YOURMOUNTDISK'``. It allows to detect others disk inside of the singularity container. $YOURMOUNTDISK corresponds to mount disk, it could be $HOME or another disk path.

.. note::
    For others snakemake arguments, please check documentation https://snakemake.readthedocs.io/en/v5.11.0/executing/cli.html#all-options


Cluster execution
-----------------

This is a typical launcher for using CulebrONT on a SLURM cluster. You have to adapt it for the configuration of your favorite one. You can use wrappers or profiles.

wrappers
~~~~~~~~

A ``slurm_wrapper.py`` script is available on CulebrONT projet to manage resources from your cluster configuration (from cluster_config.yaml file). This is the easier way to know what is running on cluster and to adapt resources for every job. Take care, this cluster_config.yaml file is becoming obsolete on latest Snakemake versions.

.. code-block::bash

   #!/bin/bash
   #SBATCH --job-name culebrONT
   #SBATCH --output slurm-%x_%j.log
   #SBATCH --error slurm-%x_%j.log

   module load system/singularity/3.3.0
   module load system/python/3.7.2

   snakemake --unlock

   # SLURM JOBS WITH USING WRAPPER
   snakemake --nolock --use-conda --use-singularity --cores -p -s Snakefile --latency-wait 60000000 --keep-going --restart-times 0 --rerun-incomplete --configfile config.yaml --cluster "python3 slurm_wrapper.py config.yaml cluster_config.yaml" --cluster-config cluster_config.yaml --cluster-status "python3 slurm_status.py"


profiles
~~~~~~~~

Optionally is possible to use Profiles in order to run CulebrONT on HPC cluster. Please follow the `recommendations found on the SnakeMake profile github <https://github.com/Snakemake-Profiles/>`_.

.. code-block::bash

   #!/bin/bash
   #SBATCH --job-name culebrONT
   #SBATCH --output slurm-%x_%j.log
   #SBATCH --error slurm-%x_%j.log

   module load system/singularity/3.3.0
   module load system/python/3.7.2

   snakemake --unlock

   # USING PROFILES
   snakemake --nolock --use-singularity --use-conda --cores -p -s Snakefile --configfile config.yaml --latency-wait 60000000 --keep-going --restart-times 0 --rerun-incomplete --cluster-config cluster_config.yaml --profile slurm-culebrONT

.. note::
    For others snakemake cluster arguments, please check documentation https://snakemake.readthedocs.io/en/stable/executing/cluster.html


In any case, this launcher can be submitted to the SLURM queue typing:

.. code-block::

    sbatch submit_culebront.sh



Output on CulebrONT
===================

The architecture of CulebrONT output is designed as follows:

.. code-block::bash

    OUTPUT_CULEBRONT_CIRCULAR/
    ├── SAMPLE-1
    │   ├── AGGREGATED_QC
    │   │   ├── DATA
    │   │   ├── MAUVE_ALIGN
    │   │   └── QUAST_RESULTS
    │   ├── ASSEMBLERS
    │   │   ├── CANU
    │   │   │   ├── ASSEMBLER
    │   │   │   ├── CORRECTION
    │   │   │   ├── FIXSTART
    │   │   │   ├── POLISHING
    │   │   │   └── QUALITY
    │   │   ├── FLYE
    │   │   │   ├── ...
    │   │   ├── MINIASM
    │   │   │   ├── ...
    │   │   ├── RAVEN
    │   │   │   ├── ...
    │   │   ├── SHASTA
    │   │   │   ├── ...
    │   │   └── SMARTDENOVO
    │   │   │   ├── ...
    │   ├── DIVERS
    │   │   └── FASTQ2FASTA
    │   ├── LOGS
    │   └── REPORT
    └── FINAL_REPORT
    ├── SAMPLE-2 ...


Report
======

CulebrONT generates a useful report containing, foreach fastq, a summary of interesting statistics. Please discover an |location_link| ... and enjoy !!

.. |location_link| raw:: html

    <a href="https://itrop.ird.fr/culebront_utilities/FINAL_REPORT/CulebrONT_report.html" target="_blank">example</a>


.. important::
    To visualise the report created by CulebrONT, transfer the folder ``FINAL_RESULTS`` on your local computer and open it on a navigator.

