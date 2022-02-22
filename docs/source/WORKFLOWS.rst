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

To create this file juste run


.. click:: culebrONT.main:create_config
    :prog: culebrONT create_config
    :show-nested:

Then edit the relevant sections of the file to customize your flavor of a workflow.

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
    "FAST5","Nanopolish uses FAST5 files for polishing and Medaka needs FAST5 files if a model training step is requested. Please give the path of the FAST5 folder in the *FAST5* DATA parameter. Inside this directory, a subdirectory with the exact same name as the corresponding FASTQ (before the *.fastq.gz*\ ) is required. For instance, if in the *FASTQ* directory we have *run1.fastq.gz* and *run2.fastq.gz*\ , CulebrONT is expecting the *run1/* and *run2/* subdirectories in the FAST5 main directory"
    "ILLUMINA","Indicate the path to the directory with *Illumina* sequence data (in fastq or fastq.gz format) to perform pilon correction and for KAT on quality. Use preferentially paired-end data. All fastq files need to be homogeneous on their extension name. Please use *run1_R1* and *run1_R2* nomenclature."
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
    If you expect your genome to include a circular replicon (e.g. with prokaryote), it is recommendated to activate CIRCULAR steps

Example:

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 10-27


3. Choose quality tools
-----------------------
With CulebrONT you can use several quality tools to check assemblies.

* If BUSCO or QUAST are used, they will run on every fasta assembly generated along the various steps of the pipeline.

* If BLOBTOOLS, ASSEMBLYTICS, FLAGSTATS and KAT are activated only the fasta assembly generated after the last sequence processing step of the pipeline will be used.

* KAT quality tool can be activate but Illumina reads are mandatory in this case. These reads can be compressed or not.

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 29-38

If several assemblers are activated, a multiple alignment of the various assemblies **for small genomes (<10-20Mbp)** can be computed with Mauve.

* If you want to improve alignment with MAUVE on circular molecules, it is recommended to activate the *Fixstart* step.
* Only activate MAUVE if you have more than one assembler per sample and more than one quality step.

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 41-43


4. Parameters for some specific tools
--------------------------------------

You can manage tools parameters on the params section on ``config.yaml`` file.

Specifically for ``Racon``:

* Racon can be launched recursively from 1 to 9 rounds.

Specifically for ``Medaka``:

* If 'MEDAKA_TRAIN_WITH_REF' is activated, Medaka launchs training using the reference found in 'DATA/REF' param. Medaka does not take into account other medaka model parameters and use the trained model instead.

* If 'MEDAKA_TRAIN_WITH_REF' is deactivated, Medaka does not launch training but uses instead the model provided in 'MEDAKA_MODEL_PATH'. Give to CulebrONT the path of medaka model *OR* just the model name in order to correct assemblies. This parameter could not be empty.

.. important::
    Medaka models can be downloaded from the medaka repository. You need to install ``git lfs`` (see documentation here https://git-lfs.github.com/) to download largest files before ``git clone https://github.com/nanoporetech/medaka.git\``.

Specifically for ``Pilon``:

* We set java memory into Singularity.culebront_tools to 8G. If you need to allocate more memory it's possible to chang this line ``sed -i "s/-Xmx1g/-Xmx8g/g" /usr/local/miniconda/miniconda3/envs/pilon/bin/pilon`` in the ``Containers/Singularity.culebront_tools.def`` recipe file before building the Singularity image.

Specifically for ``Busco``:

* If BUSCO is activated, one must provide CulebrONT the path of busco a database *OR* only the database name (See the `Busco documentation <https://busco.ezlab.org/busco_userguide.html#genome-mode-assessing-a-genome-assembly>`_).This parameter cannot be empty.

Specifically for ``Blobtools``:
* nodes and names from the ncbi taxdump database can be download from here : https://github.com/DRL/blobtools#download-ncbi-taxdump-and-create-nodesdb

Here you find standard parameters used on CulebrONT. Feel free to adapt it to your own requirements.

.. literalinclude:: ../../culebrONT/install_files/config.yaml
    :language: YAML
    :lines: 46-

.. warning::
    Please check documentation of each tool and make sure that the settings are correct!


------------------------------------------------------------------------

How to run the workflow
=======================

Before attempting to run CulebrONT please be sure you have already modified the ``config.yaml`` file as explained in :ref:`1. Providing data`.

If you installed CulebrONT on a HPC cluster with a job scheduler, you can run:


.. click:: culebrONT.main:run_cluster
    :prog: culebrONT run_cluster
    :show-nested:


------------------------------------------------------------------------


.. click:: culebrONT.main:run_local
    :prog: culebrONT run_local
    :show-nested:

------------------------------------------------------------------------

Advance run
===========

Provide more ressources
-----------------------

If cluster default resources are not sufficient, you can edit ``cluster_config.yaml`` See :ref:`2. Adapting *cluster_config.yaml*`:

.. click:: culebrONT.main:edit_cluster_config
    :prog: culebrONT edit_cluster_config
    :show-nested:

Now, take a coffee or tea, and enjoy !!!!!!

------------------------------------------------------------------------

Provide own tools_config.yaml
-----------------------------

To change the tools used on culebrONT workflow, you can run See :ref:`3. How to configure tools_path.yaml`

.. click:: culebrONT.main:edit_tools
    :prog: culebrONT edit_tools
    :show-nested:

------------------------------------------------------------------------

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

CulebrONT generates a useful report including the versions of tools used and, foreach fastq, a summary of interesting statistics and . Please discover an |location_link| ... and enjoy !!

.. note::

    Because of constraints imposed by Snakemake, we have not been able to include the version of bwa and seqtk in the report https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-5-loggin. If you want to know the versions of these tools, go check by yourself ^^.


.. |location_link| raw:: html

    <a href="https://itrop.ird.fr/culebront_utilities/FINAL_REPORT/CulebrONT_report.html" target="_blank">example</a>


.. important::

    To visualise the report created by CulebrONT, transfer the folder ``FINAL_RESULTS`` on your local computer and open it on a navigator.

