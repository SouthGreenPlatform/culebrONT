.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Requirements
============

CulebrONT requires |PythonVersions| and |graphviz|.

CulebrONT has been mostly developed to work on an HPC but a local installation is also possible.

------------------------------------------------------------------------

Install culebrONT PyPI package
===============================

First install culebrONT python package with pip.

.. code-block:: bash

   python3 -m pip install culebrONT
   culebrONT --help

Now, follow this documentation according to what you want, local or HPC mode.

------------------------------------------------------------------------

Steps for LOCAL installation
============================

Install CulebrONT in a *local* mode using ``culebrONT install_local`` command line.

.. click:: culebrONT.main:install_local
    :prog: culebrONT install_local
    :show-nested:

To create a pipeline, tools used by CulebrONT are wrapped into ``Singularity images``. These images are automatically downloaded and used to configure files needed by the pipeline.  Local mode install, without scheduler, is constrains to use singularity.


After install, optionally (but recommended) you can now check CulebrONT installation using a scaled dataset.
See the section :ref:`Check install` for details.

------------------------------------------------------------------------

Steps for HPC installation
==========================

CulebrONT uses available snakemake profiles to easier cluster installation and resources management.
Run the command `culebrONT install_cluster` to install on HPC.
We tried to make HPC installation as easy as possible but it's necessary to adapt some files according to your HPC environment.


.. click:: culebrONT.main:install_cluster
    :prog: culebrONT install_cluster
    :show-nested:

1. Adapt `profile` and `cluster_config.yaml`
---------------------------------------------

Now that CulebrONT is installed, it offers default configuration files but they can be modified. Please check and adapt these files to your system architecture.

1. Adapt the pre-formatted `f –env si`snakemake profile`` to configure cluster options.
See the section :ref:`1. Snakemake profiles` for details.

2. Adapt the :file:`cluster_config.yaml` file to manage cluster resources such as partition, memory and threads available for each job.
See the section :ref:`2. Adapting *cluster_config.yaml*` for further details.


2. Adapt `tools_path.yaml`
--------------------------

As CulebrONT uses many tools, you must install them using one of two possibilities:

1. Either through the |Singularity| containers,

2. Or using the ``module load`` mode,

.. code-block:: bash

   culebrONT install_cluster --help
   culebrONT install_cluster --scheduler slurm --env modules
   # OR
   culebrONT install_cluster --scheduler slurm --env singularity

If ``--env singularity`` argument is used, CulebrONT download previously build images containing environment need to run CulebrONT (tools and dependencies).

Adapt the file :file:``tools_path.yaml`` - in YAML (Yet Another Markup Language) - format  to indicate to CulebrONT where the tools are installed.
See the section :ref:`3. How to configure tools_path.yaml` for details.


------------------------------------------------------------------------

Check install
==============

Optionally, in order to test your install of CulebrONT pipeline, a data test called ``Data-Xoo-sub/`` is available on https://itrop.ird.fr/culebront_utilities/.

.. click:: culebrONT.main:test_install
    :prog: culebrONT test_install
    :show-nested:

This dataset is automatically downloaded by culebrONT in the ``-d`` repertory using :

.. code-block:: bash

   culebrONT test_install -d test

Launching suggested command line done by CulebrONT, in CLUSTER mode :

.. code-block:: bash

    culebrONT run_cluster --config test/data_test_config.yaml

Or in local mode :

.. code-block:: bash

    culebrONT run_local -t 8 -c test/data_test_config.yaml --singularity-args "--bind $HOME"


------------------------------------------------------------------------

Advance installation
====================


1. Snakemake profiles
---------------------

The Snakemake-profiles project is an open effort to create configuration profiles allowing to execute Snakemake in various computing environments
(job scheduling systems as Slurm, SGE, Grid middleware, or cloud computing), and available at https://github.com/Snakemake-Profiles/doc.

In order to run CulebrONT on HPC cluster, we use profiles.

Quickly, see `here <https://github.com/SouthGreenPlatform/culebrONT/blob/master/culebrONT/install_files/cluster_config_SLURM.yaml>`_ an example of the Snakemake SLURM profile we use for the French national bioinformatics infrastructure at IFB.

More info about profiles can be found here https://github.com/Snakemake-Profiles/slurm#quickstart.

Preparing the profile's *config.yaml* file
******************************************

Once, your basic profile is created. To finalize it, modify as necessary the ``culebrONT/culebrONT/default_profile/config.yaml`` to customize Snakemake parameters that will be used internally by CulebrONT:

.. code-block:: ini

    restart-times: 0
    jobscript: "slurm-jobscript.sh"
    cluster: "slurm-submit.py"
    cluster-status: "slurm-status.py"
    max-jobs-per-second: 1
    max-status-checks-per-second: 10
    local-cores: 1
    jobs: 200                   # edit to limit number jobs submit
    latency-wait: 60000000
    use-envmodules: true        # Adapt True/False only active one
    use-singularity: false      # if False, please install all R package on tools_config.yaml ENVMODULE/R
    rerun-incomplete: true
    printshellcmds: true


2. Adapting *cluster_config.yaml*
----------------------------------

In the ``cluster_config.yaml`` file, you can manage HPC resources, choosing partition, memory and threads to be used by default or specifically for each rule/tool depending on your HPC Job Scheduler (see `there <https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration-deprecated>`_). This file generally belongs to a Snakemake profile (see above).

.. warning::
    If more memory or threads are requested, please adapt the content of this file before running on your cluster.

A list of CulebrONT rules names can be found in the section :ref:`Threading rules inside culebrONT`

.. warning::
    For some rules in the *cluster_config.yaml* as `rule_graph` or `run_get_versions` we use by default wildcards, please don't remove it.


3. How to configure tools_path.yaml
-----------------------------------

In the ``tools_path`` file, you can find two sections: SINGULARITY and ENVMODULES. In order to fill it correctly, you have 2 options:

1. Use only SINGULARITY containers. In this case, fill only this section. Put the path to the build singularity images you want to use.
Absolute paths are strongly recommended but not mandatory. See the section :ref:`'How to build singularity images'<How to build singularity images>`  for further details.

.. literalinclude:: ../../culebrONT/install_files/tools_path.yaml
    :language: YAML
    :lines: 6-8

.. warning::
    For SINGULARITY containers to be actually used, one needs to make sure that the *--use-singularity* flag is included in the snakemake command line.

2. Use only ENVMODULES. In this case, fill this section with modules available on your cluster (here is an example):

.. literalinclude:: ../../culebrONT/install_files/tools_path.yaml
    :language: YAML
    :lines: 10-18

CulebrONT needs a wide set of R modules for reporting, if you use ENVMODULE R. Just have a look at dependencies in the ``Containers/Singularity.report.def`` file.
Yes, plenty of packages!! That's why we provide build singularity containers ready to use and recommend them for the R part.

.. note::

    TIP !! We provide a singularity container for R packages (Singularity.report.def), you can use this one to create a module environment.

.. warning::
    Make sure to specify the *--use-envmodules* flag in the snakemake command line for ENVMODULE to be implemented.
    More details can be found here: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules


------------------------------------------------------------------------

And more ...
-------------

How to build singularity images
*******************************

You can build your own image using the available *.def* recipes from the ``culebrONT/culebrONT/containers/`` directory.

.. warning::
    Be careful, you need root access to build singularity images

.. code-block:: bash

    cd culebrONT/culebrONT/containers/
    sudo make build

Threading rules inside culebrONT
********************************

Please find here the rules names found in CulebrONT code.
It is recommended to set threads using the snakemake command when running on a single machine or in the cluster configuration file to manage cluster resources via the job scheduler.
This would save the user a painful exploration of the snakefiles code of CulebrONT.

.. code-block:: python

    run_flye
    run_canu
    run_minimap_for_miniasm
    run_miniasm
    run_minipolish
    run_raven
    convert_fastq_to_fasta
    run_smartdenovo
    run_shasta
    run_circlator
    tag_circular
    tag_circular_to_minipolish
    rotate_circular
    run_fixstart
    run_makerange
    run_nanopolish_index
    preparing_ref_to_nanopolish
    run_nanopolish_variants
    run_nanopolish_merge
    index_fasta_to_correction
    run_minialign_to_medaka
    run_medaka_train
    run_medaka_consensus
    run_medaka_merge
    run_pilon_first_round
    run_pilon
    run_racon
    preparing_fasta_to_quality
    run_quast
    run_busco
    run_diamond
    run_minimap2
    run_blobtools
    run_mummer
    run_assemblytics
    combined_fastq
    run_KAT
    run_mauve
    run_bwa_mem2
    run_flagstat
    final
    rule_graph
    run_report_snakemake
    run_flagstats_stats
    run_busco_stats
    run_racon_version
    run_busco_version
    run_benchmark_time
    run_get_versions
    stats_assembly
    run_report



.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7%2B-blue
   :target: https://www.python.org/downloads
   :alt: Python 3.7+

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake 5.10.0+

.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/
   :alt: Singularity 3.10.0+

.. |graphviz| image:: https://img.shields.io/badge/graphviz-%3E%3D2.40.1-green
   :target: https://graphviz.org/
   :alt: graphviz 2.40.1+
