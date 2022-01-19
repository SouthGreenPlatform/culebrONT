.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Requirements
============

CulebrONT requires |PythonVersions|, |SnakemakeVersions| and |graphviz|.

CulebrONT has been mostly developed to work on an HPC but a local installation is also possible.

------------------------------------------------------------------------

Steps for LOCAL installation
============================

CulebrONT and dependencies and also every tool used to create a pipeline are available through a ``Singularity images``.
To install CulebrONT on *local* mode you must use singularity:

1. Install culebrONT package
----------------------------

First install culebrONT python package with pip.

.. code-block:: bash

   python3 -m pip install culebrONT
   culebrONT --help

Then run the commande line to install on LOCAL

.. code-block:: bash

   culebrONT install_local --help
   culebrONT install_local

The script automatically download singularity images required and configure culebrONT.

Test CulebrONT install (Optional but recommended) using an available dataset.
See the section :ref:`Check install` for details.

------------------------------------------------------------------------

Steps for HPC installation
==========================

As CulebrONT uses many tools, you must install them through two possibilities:

1. Either through the |Singularity| containers,

2. Or using the ``module load`` mode,

Let's check steps for **HPC installation** :

First install culebrONT python package with pip.

.. code-block:: bash

   python3 -m pip install culebrONT
   culebrONT --help

Then run the commande line to install on HPC

.. code-block:: bash

   culebrONT install_cluster --help
   culebrONT install_cluster --scheduler slurm --env modules --bash_completion --create_envmodule --modules_dir /path/to/dir/culebrONT/

the script uses the snakemake profiles to build the installation profile for culebrONT.
if --env is singularity, culebrONT download images.
Then, the script proposes to modify the following files to adapt to your system achitecture


1. Adapt the :file:`cluster_config.yaml` file to manage cluster resources such as partition, memory and threads available for each job.
See the section :ref:`1. Preparing *cluster_config.yaml*` for further details.

2. Create a *snakemake profile* to configure cluster options.
See the section :ref:`2. Snakemake profiles` for details.

3. Adapt the file :file:`tools_path.yaml` - in YAML (Yet Another Markup Language) - format  to indicate to CulebrONT where the tools are installed.
See the section :ref:`3. How to configure tools_path.yaml` for details.

4. Test CulebrONT install (Optional but recommended) using an available dataset.
See the section :ref:`Check install` for details.



1. Preparing *cluster_config.yaml*
----------------------------------

In the ``cluster_config.yaml`` file, you can add partition, memory and threads to be used by default or specifically for each rule/tool.


.. warning::
    If more memory or threads are requested, please adapt the content of this file before running on your cluster.

Here is a example of the configuration file we used on our `i-Trop HPC <../../cluster_config.yaml>`_ .

A list of rules names can be found in the section :ref:`Threading rules inside culebrONT`

.. warning::
    Please give to *cluster_config.yaml* specific parameters to rules get_versions and rule_graph without using wildcards into log files.


2. Snakemake profiles
---------------------

The Snakemake-profiles project is an open effort to create configuration profiles allowing to execute Snakemake in various computing environments
(job scheduling systems as Slurm, SGE, Grid middleware, or cloud computing), and available at https://github.com/Snakemake-Profiles/doc.

In order to run CulebrONT on HPC cluster, we use profiles.

Quickly, here is an example of the Snakemake SLURM profile we use for the French national bioinformatics infrastructure at IFB.
We followed the documentation found here https://github.com/Snakemake-Profiles/slurm#quickstart.

Now, your basic profile is created. To finalize it, change the ``CulebrONT_pipeline/profiles/CulebrONT/config.yaml`` to :

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


.. note::
   You can find the generated files when profiles are created following this documentation. IFB profile example can be found on the `̀ gift_files/CulebrONT_SLURM/`` repertory on our github repository.

3. How to configure tools_path.yaml
-----------------------------------

In the ``|tools_path|`` file, you can find two sections: SINGULARITY and ENVMODULES. In order to fill it correctly, you have 2 options:

1. Use only SINGULARITY containers. In this case, fill only this section. Put the path to the builded singularity images.
Absolute paths are strongly recommended but not mandatory. See the section :ref:`'How to build singularity images'<How to build singularity images>`  for further details.

.. literalinclude:: ../../culebrONT/install_files/tools_path.yaml
    :language: YAML
    :lines: 6-8

.. warning::
    The use of SINGULARITY is constraint with the use the *--use-singularity* parameter in the snakemake command line.

2. Use only ENVMODULES. In this case, fill this section with modules available on your cluster (here is an example):

.. literalinclude:: ../../culebrONT/install_files/tools_path.yaml
    :language: YAML
    :lines: 10-18

CulebrONT needs a wide set of R modules for reporting, if you use ENVMODULE R. Just have a look at dependencies in the ``Containers/Singularity.report.def`` file.
Yes, plenty of packages!! That's why we provide build singularity containers ready to use and recommend them for the R part.

.. note::

    TIP !! We provide a singularity container for R packages (Singularity.report.def), you can use this one to create a module environment.

.. warning::
    The use of ENVMODULE is constraint with the use the *--use-envmodules* parameter in the snakemake command line.
    More details can be found here: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules

------------------------------------------------------------------------

Check install
=============

Optionally, in order to test your install of CulebrONT pipeline, a data test called ``Data-Xoo-sub/`` is available on https://itrop.ird.fr/culebront_utilities/.
Feel free to download it using ``wget`` and put it on CulebrONT directory.

.. code-block:: bash

    cd CulebrONT_pipeline
    wget --no-check-certificat -rm -nH --cut-dirs=1 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/Data-Xoo-sub/


Now, it is time to prepare the configuration file ``config.yaml`` file to indicate to CulebrONT which kind of pipeline you want to create and to run on your data ! :ref:`How to create a workflow`.

Finally run CulebrONT!!

.. code-block:: bash

    culebrONT run_cluster --config config.yaml --dry-run
    culebrONT run_local --config config.yaml --threads 6 --dry-run

------------------------------------------------------------------------

Advance installation
====================

How to build singularity images
-------------------------------

You can build your own image using the available *.def* recipes from the ``culebrONT/culebrONT/containers/`` directory.

.. warning::
    Be careful, you need root rights to build singularity images

.. code-block:: bash

    cd culebrONT/culebrONT/containers/
    sudo make build

------------------------------------------------------------------------

Threading rules inside culebrONT
--------------------------------

Please find here the rules names found in CulebrONT code.
It could be useful to set threads in local running using the snakemake command or in the cluster configuration to manage cluster resources.
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
