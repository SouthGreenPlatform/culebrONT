.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

How to create a workflow
========================

CulebrONT allow you to build a workflow using a simple configuration ``config.yaml`` file. In this file :

* First, provide data paths
* Second, activate tools from assembly to correction.
* Third, activate tools from quality checking of assemblies.
* And last, manage parameters tools.

To create file juste run

.. code-block:: bash

   culebrONT create_config --help
   culebrONT create_config -configyaml P/ath/to/file

Then edit section on file according to your workflow.

1. Providing data
------------------

First, indicate the data path in the configuration ``config.yaml`` file:

.. code-block:: YAML

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

    "FASTQ", "Every FASTQ file should contain the whole set of reads to be assembled. Each fastq file will be assembled independently."
    "REF","Only one REFERENCE genome file will be used by CulebrONT. This REFERENCE will be used for quality steps (ASSEMBLYTICS, QUAST and MAUVE)"
    "GENOME_SIZE", "Estimated genome size of the assembly can be done on mega (Mb), giga(Gb) or kilobases (Kb). This size is used on some assemblers (CANU) and also on QUAST quality step"
    "FAST5","Nanopolish needs FAST5 files to training steps. Please give the path of FAST5 folder in the *FAST5* DATA parameter. Inside this directory, a subdirectory with the exact same name as the corresponding FASTQ (before the *.fastq.gz*\ ) is requested. For instance, if in the *FASTQ* directory we have *run1.fastq.gz* and *run2.fastq.gz*\ , CulebrONT is expecting the *run1/* and *run2/* subdirectories in the FAST5 main directory"
    "ILLUMINA","Indicate the path to the directory with *Illumina* sequence data (in fastq or fastq.gz format) to perform pilon correction and KAT on quality. Use preferentially paired-end data. All fastq files need to be homogeneous on their extension name. Please use *sample_R1* and *sample_R2* nomenclature."
    "OUTPUT","output *path* directory"

.. warning::

    For FASTQ, naming convention accepted by CulebrONT is *NAME.fastq.gz* or *NAME.fq.gz* or *NAME.fastq* or *NAME.fq*. Preferentially use short names and avoid special characters because report can fail. Avoid to use the long name given directly by sequencer.

    All fastq files have to be homogeneous on their extension and can be compressed or not.

    Reference fasta file need a fasta or fa extension uncompressed.

2. Choose assemblers, polisher and correctors
---------------------------------------------

Activate/deactivate assemblers, polishers and correctors as you wish.
Feel free to activate only assembly, assembly+polishing or assembly+polishing+correction.

.. note::
    If you are working on prokaryote, is recommendated to activate CIRCULAR steps

Example:

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 10-27


3. Choose quality tools
-----------------------
With CulebrONT you can use several quality tools to check assemblies.

* If BUSCO or QUAST are used, every fasta generated on the pipeline will be used with them.

* If BLOBTOOLS, ASSEMBLYTICS, FLAGSTATS and KAT are activated only the last draft generated on the pipeline will be used.

* KAT quality tool can be activate but Illumina reads are mandatory in this case. These reads can be compressed or not.

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 29-38

Alignment of various assemblies **for small genomes (<10-20Mbp)** is also possible by using Mauve.

* If you want to improve alignment with MAUVE on circular molecules is recommended to activate *Fixstart* step.
* Only activate MAUVE if you have more than one sample and more than one quality step.

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 41-43


4. Parameters for some specific tools
--------------------------------------

You can manage tools parameters on the params section on ``config.yaml`` file.

Specifically to ``Racon``:

* Racon can be launch recursively from 1 to 9 rounds. 2 or 3 are recommended.

Specifically to ``Medaka``:

* If 'MEDAKA_TRAIN_WITH_REF' is activated, Medaka launchs training using the reference found in 'DATA/REF' param. Medaka does not take into account other medaka model parameters.

* If 'MEDAKA_TRAIN_WITH_REF' is deactivated, Medaka does not launch training but uses instead the model provided in 'MEDAKA_MODEL_PATH'. Give to CulebrONT path of medaka model *OR* only the model name in order to correct assemblies. This parameter could not be empty.

.. important::
    Medaka models can be downloaded from the medaka repository. You need to install ``git lfs`` (see documentation here https://git-lfs.github.com/) to download largest files before ``git clone https://github.com/nanoporetech/medaka.git\``.

Specifically to ``Pilon``:

* We set java memory into Singularity.culebront_tools to 8G. If you need to allocate more memory it's possible changing this line ``sed -i "s/-Xmx1g/-Xmx8g/g" /usr/local/miniconda/miniconda3/envs/pilon/bin/pilon`` into ``Containers/Singularity.culebront_tools.def`` before singularity building images.

Specifically to ``Busco``:

* If BUSCO is activate, give to CulebrONT the path of busco database *OR* only the database name.This parameter could not be empty.

Specifically to ``Blobtools``:
* nodes and names from ncbi taxdump database can be download from here : https://github.com/DRL/blobtools#download-ncbi-taxdump-and-create-nodesdb

Here you find standard parameters used on CulebrONT. Feel free to adapt it to your requires.

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 46-

.. warning::
    Please check documentation of each tool and make sure that the settings are correct!

.. ############################################################

How to run the workflow
=======================

Before run CulebrONT please be sure you have already modified the ``config.yaml`` file as was explained on :ref:`1. Providing data`

.. code-block:: python

    culebrONT run_cluster --help
    culebrONT run_cluster -c config.yaml --dry-run


You can optionally also pass to snakemake more options by using non option parameter (check it https://snakemake.readthedocs.io/en/stable/executing/cli.html).

.. code-block:: bash

    # in LOCAL using maximum 8 threads
    culebrONT run_local -c config.yaml --cores 8 --dry-run

    # in LOCAL using 6 threads for Canu assembly from the total 8 threads
    culebrONT run_local -c config.yaml ---cores 8 --set-threads run_canu=6


Expert mode
===========

Provide more ressources
-----------------------

If cluster default resources are not sufficient, you can overwrite ``cluster_config.yaml`` file used by the profile doing:

.. code-block:: bash

    # in HPC overwriting cluster_config.yaml given by user
    culebrONT create_cluster_config --clusterconfig own_cluster_params.yaml
    culebrONT run_cluster --config config.yaml --clusterconfig own_cluster_params.yaml

Now, take a coffee or tea, and enjoy !!!!!!

Provide own tools_config.yaml
-----------------------------

To change the tools used on culebrONT workflow, you can run

.. code-block:: bash

    culebrONT edit_tools


Output on CulebrONT
===================

The architecture of CulebrONT output is designed as follows:

.. code-block:: bash

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

CulebrONT generates a useful report containing, foreach fastq, a summary of interesting statistics and versions of tools used. Please discover an |location_link| ... and enjoy !!

.. note::
    Because of constraints imposed by Snakemake we have not been able to recover the version of bwa and seqtk on report https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-5-loggin. If you want to know the versions of these tools, go check by yourself ^^.


.. |location_link| raw:: html

    <a href="https://itrop.ird.fr/culebront_utilities/FINAL_REPORT/CulebrONT_report.html" target="_blank">example</a>


.. important::
    To visualise the report created by CulebrONT, transfer the folder ``FINAL_RESULTS`` on your local computer and open it on a navigator.

