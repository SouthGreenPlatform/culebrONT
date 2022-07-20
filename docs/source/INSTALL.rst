.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Requirements
============

CulebrONT requires |PythonVersions|, |graphviz| and |RVersions|.

CulebrONT is developed to work mostly on an HPC distributed cluster but a local, single machine, installation is also possible.

------------------------------------------------------------------------

Install CulebrONT PyPI package
===============================

First, install the CulebrONT python package with pip.

.. code-block:: bash

   python3 -m pip install culebrONT
   culebrONT --help

Now, follow this documentation according to what you want, local or HPC mode.

------------------------------------------------------------------------

Steps for LOCAL installation
============================

Install CulebrONT in a *local* (single machine) mode using ``culebrONT install_local`` command line.

.. click:: culebrONT.main:install_local
   :prog: culebrONT install_local
   :show-nested:

To create a pipeline, tools used by CulebrONT are wrapped into ``Singularity images``. These images are automatically downloaded and used by the configuration files of the pipeline. Local mode install, without scheduler, is constrains to use these Singularity images.

.. warning::
   Singularity images are downloaded in the location of the package CulebrONT. Be careful these images need at approximatly 4G of free space. If installed with Pypi with the flag --user (without root), the package is installed in your HOME.


Optionally (but recommended), after installing in local, you can check the CulebrONT installation using a dataset scaled for single machine.
See the section :ref:`Check install` for details.

------------------------------------------------------------------------

Steps for HPC distributed cluster installation
==============================================

CulebrONT uses any available snakemake profiles to ease cluster installation and resources management.
Run the command `culebrONT install_cluster` to install on a HPC cluster.
We tried to make cluster installation as easy as possible, but it is somehow necessary to adapt a few files according to your cluster environment.


.. click:: culebrONT.main:install_cluster
   :prog: culebrONT install_cluster
   :show-nested:

1. Adapt `profile` and `cluster_config.yaml`
---------------------------------------------
f
Now that CulebrONT is installed, it proposes default configuration files, but they can be modified. Please check and adapt these files to your own system architecture.

1. Adapt the pre-formatted `f –env si`snakemake profile`` to configure your cluster options.
See the section :ref:`1. Snakemake profiles` for details.

2. Adapt the :file:`cluster_config.yaml` file to manage cluster resources such as partition, memory and threads available for each job.
See the section :ref:`2. Adapting *cluster_config.yaml*` for further details.


2. Adapt `tools_path.yaml`
--------------------------

As CulebrONT uses many tools, you must install them using one of the two following possibilities:

1. Either through the |Singularity| containers,

2. Or using the ``module load`` mode,

.. code-block:: bash

   culebrONT install_cluster --help
   culebrONT install_cluster --scheduler slurm --env modules
   # OR
   culebrONT install_cluster --scheduler slurm --env singularity

If ``--env singularity`` argument is specified, CulebrONT will download previously build Singularity images, containing the complete environment need to run CulebrONT (tools and dependencies).

Adapt the file :file:``tools_path.yaml`` - in YAML (Yet Another Markup Language) - format  to indicate CulebrONT where the different tools are installed on your cluster.
See the section :ref:`3. How to configure tools_path.yaml` for details.


------------------------------------------------------------------------

Check install
==============

In order to test your install of CulebrONT, a data test called ``Data-Xoo-sub/`` is available at https://itrop.ird.fr/culebront_utilities/.

.. click:: culebrONT.main:test_install
   :prog: culebrONT test_install
   :show-nested:

This dataset will be automatically downloaded by CulebrONT in the ``-d`` repertory using :

.. code-block:: bash

   culebrONT test_install -d test

Launching the (suggested, to be adapted) command line in CLUSTER mode will perform the tests:

.. code-block:: bash

   culebrONT run_cluster --config test/data_test_config.yaml

In local mode, type :

.. code-block:: bash

   culebrONT run_local -t 8 -c test/data_test_config.yaml --singularity-args "--bind $HOME"


------------------------------------------------------------------------

Advance installation
====================


1. Snakemake profiles
---------------------

The Snakemake-profiles project is an open effort to create configuration profiles allowing to execute Snakemake in various computing environments
(job scheduling systems as Slurm, SGE, Grid middleware, or cloud computing), and available at https://github.com/Snakemake-Profiles/doc.

In order to run CulebrONT on HPC cluster, we take advantages of profiles.

Quickly, see `here <https://github.com/SouthGreenPlatform/culebrONT/blob/master/culebrONT/install_files/cluster_config_SLURM.yaml>`_ an example of the Snakemake SLURM profile we used for the French national bioinformatics infrastructure at IFB.

More info about profiles can be found here https://github.com/Snakemake-Profiles/slurm#quickstart.

Preparing the profile's *config.yaml* file
******************************************

Once your basic profile is created, to finalize it, modify as necessary the ``culebrONT/culebrONT/default_profile/config.yaml`` to customize Snakemake parameters that will be used internally by CulebrONT:

.. code-block:: ini

   restart-times: 0
   jobscript: "slurm-jobscript.sh"
   cluster: "slurm-submit.py"
   cluster-status: "slurm-status.py"
   max-jobs-per-second: 1
   max-status-checks-per-second: 10
   local-cores: 1
   jobs: 200                   # edit to limit the number of jobs submitted in parallel
   latency-wait: 60000000
   use-envmodules: true        # adapt True/False for env of singularuty, but only active one possibility !
   use-singularity: false      # if False, please install all R packages listed in tools_config.yaml ENVMODULE/R
   rerun-incomplete: true
   printshellcmds: true


2. Adapting *cluster_config.yaml*
----------------------------------

In the ``cluster_config.yaml`` file, you can manage HPC resources, choosing partition, memory and threads to be used by default,
or specifically, for each rule/tool depending on your HPC Job Scheduler (see `there <https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration-deprecated>`_). This file generally belongs to a Snakemake profile (see above).

.. warning::
   If more memory or threads are requested, please adapt the content
   of this file before running on your cluster.


A list of CulebrONT rules names can be found in the section :ref:`Threading rules inside culebrONT`


.. warning::
   For some rules in the *cluster_config.yaml* as `rule_graph` or `run_get_versions`,
   we use by default wildcards, please don't remove it.


3. How to configure tools_path.yaml
-----------------------------------

.. note::
    About versions of tools, the user can choose themself what version of tools to use with modules or with singularity.
    HOWEVER, the pipeline was validated with specific versions (check the `singularity def <https://github.com/SouthGreenPlatform/culebrONT/blob/master/culebrONT/containers/Singularity.culebront_tools.def>`_) so it may leads to error due to parameter changes.
    :ref:`Assembly`


In the ``tools_path`` file, you can find two sections: SINGULARITY and ENVMODULES. In order to fill it correctly, you have 2 options:

1. Use only SINGULARITY containers: in this case, fill only this section. Put the path to the built Singularity images you want to use.
Absolute paths are strongly recommended. See the section :ref:`'How to build singularity images'<How to build singularity images>`  for further details.

.. literalinclude:: ../../culebrONT/install_files/tools_path.yaml
   :language: YAML
   :lines: 6-8

.. warning::
   To ensure SINGULARITY containers to be really used, one needs to make
   sure that the *--use-singularity* flag is included in the snakemake command line.


2. Use only ENVMODULES: in this case, fill this section with the modules available on your cluster (here is an example):

.. literalinclude:: ../../culebrONT/install_files/tools_path.yaml
   :language: YAML
   :lines: 10-18

CulebrONT needs a wide set of R modules for reporting. If you use ENVMODULE R, just have a look at dependencies in the ``Containers/Singularity.report.def`` file.
Yes, plenty of packages!! That's why we provide build Singularity containers ready to use and recommend them for the R part.

.. note::
   TIP !! We provide a Singularity container for all R packages (Singularity.report.def),
   thus you can use this one to create a dedicated module environment.


.. warning::
   Make sure to specify the *--use-envmodules* flag in the snakemake command
   line for ENVMODULE to be implemented.
   More details can be found here:
   https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules


------------------------------------------------------------------------

And more ...
-------------

How to build Singularity images
*******************************

You can build your own image using the available *.def* recipes from the ``culebrONT/culebrONT/containers/`` directory.

.. warning::
   Be careful, you need root access to build Singularity images

.. code-block:: bash

   cd culebrONT/culebrONT/containers/
   sudo make build

Threading rules inside CulebrONT
********************************

Please find here the rules names found in CulebrONT code.
It is recommended to set threads using the snakemake command when running on a single machine,
or in a cluster configuration file to manage cluster resources through the job scheduler.
This would save users a painful exploration of the snakefiles of CulebrONT.

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

.. |RVersions| image:: https://img.shields.io/badge/R-%3E%3D4.0-red
   :target: https://cran.r-project.org/
   :alt: R 4.0+

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake 5.10.0+

.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/
   :alt: Singularity 3.10.0+

.. |graphviz| image:: https://img.shields.io/badge/graphviz-%3E%3D2.40.1-green
   :target: https://graphviz.org/
   :alt: graphviz 2.40.1+
