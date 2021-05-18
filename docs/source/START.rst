How to install CulebrONT ?
===========================

CulebrONT uses |PythonVersions| |SnakemakeVersions|  and |Singularity| as mandatory requires. In our github repository you can find available requirements (conda and singularity containers).

In your local computer OR on HPC cluster, first of all clone our repository or download last version from https://github.com/SouthGreenPlatform/CulebrONT_pipeline

.. code-block:: bash

   git clone https://github.com/SouthGreenPlatform/CulebrONT_pipeline.git
   cd CulebrONT_pipeline

In order to build a pipeline, you have to adapt two files in YAML (Yet Another Markup Language) format  ``tools_path.yaml`` and then ``config.yaml``.  The snakemake pipeline will then execute the set of rules executed into their own virtual environment.

For installation, you need to first modify ``tools_path.yaml`` file. You need to adapt path to singularity images. CulebrONT integrates many tools. In order to build a workflow, CulebrONT needs to know where tools are installed. You can give only tools you will use for your analysis. To avoid you to install them, we provide singularity images (one of them containing a lot of conda environments). We will explain all here!


How to build singularity images?
--------------------------------

You have three options to obtain singularity build images. Choose one of them!

1. Give path to **singularity hub** images from our `CulebrONT repository <https://singularity-hub.org/collections/4442>`_. If you use singularity hub repository build singularity images will be download from the cloud :

.. code-block:: YAML

    SINGULARITY:
        REPORT : 'shub://SouthGreenPlatform/CulebrONT_pipeline:report.def'
        SHASTA : 'shub://SouthGreenPlatform/CulebrONT_pipeline:shasta-0.7.0.def'
        WEESAM : 'shub://SouthGreenPlatform/CulebrONT_pipeline:weesam-latest.def'
        ASSEMBLYTICS : 'shub://SouthGreenPlatform/CulebrONT_pipeline:assemblytics-1.2.def'
        MEDAKA : 'shub://SouthGreenPlatform/CulebrONT_pipeline:medaka-gpu-1.2.def'
        BLOBTOOLS : 'shub://SouthGreenPlatform/CulebrONT_pipeline:blobtools-1.1.1.def'
        KAT: 'shub://SouthGreenPlatform/CulebrONT_pipeline:kat-latest.sif'
        BUSCO: 'shub://SouthGreenPlatform/CulebrONT_pipeline:busco-4.1.4.sif'

2. You can build available *.def* recipes available on the *CulebrONT_pipeline/Containers* repertory. Feel free to build them on your own computer (or cluster); be careful, you need root rights to do it. Use :

.. code-block:: bash

    cd Containers/
    sudo make build

Now give to CulebrONT path from build images on ``tools_path.yaml`` such as :

.. literalinclude:: ../../tools_path.yaml
    :language: YAML
    :lines: 2-11

3. Download build available singularity images from i-Trop server https://itrop.ird.fr/culebront_utilities/. Use

.. code-block:: bash

    cd CulebrONT_pipeline/Containers
    wget -rm -nH --cut-dirs=2 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/singularity_build/


.. NOTE::

    If you don't want to use some tool, please leave the path empty. Like bellow for SHASTA and ASSEMBLYTICS and don't remove the line of the tool PLEASE.

    .. code-block:: YAML

        SINGULARITY:
            REPORT : './Containers/Singularity.report.sif'
            SHASTA : ''
            WEESAM : './Containers/tools/Singularity.weesam-latest.sif'
            ASSEMBLYTICS : ''
            MEDAKA : './Containers/Singularity.medaka-gpu-1.2.sif'
            BLOBTOOLS : './Containers/Singularity.bloobtools-1.1.1.sif'
            KAT: './Containers/Singularity.kat-latest.sif'
            CONDA: './Containers/Singularity.conda.sif'


.. NOTE::

   A series of Conda environments are automatically build in the *Singularity.conda.sif* container available on itrop server. This makes installation easier. Check the Singularity.conda.def recipe to check software versions included.


HPC specifications
==================

To install CulebrONT on a global way, we recommend to change path of required dependencies on *tools_path.yaml* file and and adapt the *cluster_config.yaml* file. We recommend use singularity images as explained before.

Module load environments are usually used on clusters. If you prefer to use them, you have to fill in the *ENVMODULE* part of the cluster and don't forget to use the *--use-envmodules* parameter in the snakemake command line. More details here https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules

    .. code-block:: YAML

        ENVMODULE:
            R : "bioinfo/R/4.0.2"
            WEESAM : "bioinfo/weeSAM/1.0"
            QUAST : "bioinfo/quast/5.0.2"
            MAUVE : "bioinfo/mauve/2.4.0"
            SHASTA : "bioinfo/shasta/0.1.0"
            ASSEMBLYTICS : "bioinfo/Assemblytics/1.0"
            MEDAKA : "bioinfo/medaka-gpu/0.10.0"
            ...


Preparing *cluster_config.yaml*
-------------------------------

On ``cluster_config.yaml`` , you can add partition, memory and threads to be used by default for each rule. If more memory or threads are requested, please adapt the content of this file before running on a cluster for every rule.

.. NOTE::

   you can choose *--mem* or *--mem-per-cpu* as you see fit.

.. warning::
    please adapt the content of this file before running on a cluster for every rule !!

Here is a example of the configuration file we used on the i-Trop HPC.

.. code-block:: YAML

   __default__:
       cpus-per-task : 4
       ntasks : 1
       mem-per-cpu : '2G'
       partition : "normal"
       output : 'logs/stdout/{rule}/{wildcards}'
       error : 'logs/error/{rule}/{wildcards}'

   run_nanopolish :
       cpus-per-task : 12
       mem : '10G'
       partition : "long"

   run_canu:
       cpus-per-task : 8
       mem-per-cpu : '8G'
       partition : "long"


Available data test
===================

Optionally, in order to test install of CulebrONT pipeline, a data test ``Data-Xoo-sub/`` is available on https://itrop.ird.fr/culebront_utilities/. Feel free to download it using ``wget`` and put it on CulebrONT repertory.

.. code-block:: bash

    cd CulebrONT_pipeline
    wget -rm -nH --cut-dirs=1 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/Data-Xoo-sub/

Now, it is time to prepare configuration file ``config.yaml`` file to say to CulebrONT what kind of pipeline you want to create and use it with your data !


.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7%2B-blue
   :target: https://www.python.org/downloads
   :alt: Python 3.7+

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake 5.10.0+

.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/
   :alt: Singularity 3.10.0+

.. |Conda| image:: https://img.shields.io/badge/conda-4.8.5%20-green
   :target: https://docs.conda.io/projects/conda/en/latest/index.html
   :alt: Conda 4.8.20+
