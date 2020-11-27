How to install CulebrONT ?
===========================

CulebrONT uses |PythonVersions| |SnakemakeVersions| |Conda| and |Singularity| as mandatory requires. In our github repository you can find available conda and singularity containers.

In your local computer OR on HPC cluster, first of all clone our repository or download last version from https://github.com/SouthGreenPlatform/CulebrONT_pipeline

.. code-block:: bash

   git clone https://github.com/SouthGreenPlatform/CulebrONT_pipeline.git
   cd CulebrONT_pipeline

For installation, you need modify ``tools_path.yaml`` file. You need to adapt path to singularity images and also conda environments. We will explain all here!

How to build singularity images?
--------------------------------

You have three options to obtain singularity build images:

1. Give path to **singularity hub** images from our `CulebrONT repository <https://singularity-hub.org/collections/4442>`_. If you use singularity hub repository build singularity images will be download from the cloud :

.. code-block:: bash

    SINGULARITY:
        REPORT : 'shub://SouthGreenPlatform/CulebrONT_pipeline:report.sif'
        SHASTA : 'shub://SouthGreenPlatform/CulebrONT_pipeline:shasta-0.5.1.sif'
        WEESAM : 'shub://SouthGreenPlatform/CulebrONT_pipeline:weesam.sif'
        ASSEMBLYTICS : 'shub://SouthGreenPlatform/CulebrONT_pipeline:assemblytics.sif'
        MEDAKA : 'shub://SouthGreenPlatform/CulebrONT_pipeline:medaka-gpu-1.2.sif'
        BLOBTOOLS : 'shub://SouthGreenPlatform/CulebrONT_pipeline:bloobtools-v1.1.1.simg'
        KAT: 'path/to/Containers/tools/KAT.simg'

2. You can build available *.def* recipes available on the *CulebrONT_pipeline/Containers* repertory. Feel free to build them on your own computer (or cluster); be careful, you need root rights to do it. Use :

.. code-block:: bash

    cd Containers/
    sudo make build

Now give to CulebrONT path from build images on ``tools_path.yaml`` such as :

.. code-block:: bash

    SINGULARITY:
        REPORT : 'path/to/Containers/Singularity.report.simg'
        SHASTA : 'path/to/Containers/Singularity.shasta-0.5.1.simg'
        WEESAM : 'path/to/Containers/tools/Singularity.weesam.simg'
        ASSEMBLYTICS : 'path/to/Containers/Singularity.assemblytics-1.2.simg'
        MEDAKA : 'path/to/Containers/Singularity.medaka-gpu-1.2.simg'
        BLOBTOOLS : 'path/to/Containers/Singularity.bloobtools-v1.1.1.simg'
        KAT: 'path/to/Containers/tools/Singularity.KAT.simg'


3. Download build available singularity images from i-Trop server https://itrop.ird.fr/culebront_utilities/. Use

.. code-block:: bash

    cd CulebrONT_pipeline/Containers
    wget -rm -nH --cut-dirs=2 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/singularity_build/


How to build use Conda environments?
------------------------------------

A series of Conda environment are available on our repository for each tool used. These environments will be automatically build the first time that you launch CulebrONT. Conda are available on the *CulebrONT_pipeline/envs/* folder. Give to CulebrONT conda path required on *tools_path.yaml* file such as :

.. code-block:: bash

    CONDA:
        FLYE : './envs/flye.yaml'
        CANU : './envs/canu.yaml'
        MINIPOLISH : './envs/minipolish.yaml'
        RAVEN : './envs/raven.yaml'
        SMARTDENOVO : './envs/smartdenovo.yaml'
        CIRCULATOR : './envs/circlator.yaml'
        R : './envs/R_for_culebront_cenv.yaml'
        QUAST : './envs/quast.yaml'
        BUSCO : './envs/busco.yaml'
        DIAMOND : './envs/diamond.yaml'
        MUMMER : './envs/mummer.yaml'
        MAUVE : './envs/mauve.yaml'
        MINIASM_MINIMAP2 : './envs/miniasm_minimap2.yaml'
        MINIMAP2_SAMTOOLS : './envs/minimap2_samtools.yaml'
        RACON_MINIMAP2 : './envs/racon_minimap2.yaml'
        NANOPOLISH_MINIMAP2_SAMTOOLS_SEQTK : './envs/nanopolish_minimap2_samtools_seqtk.yaml'

.. DANGER::
    conda enviroments are compiled by Snakemake in each output analysis folder. To avoid this, please use ``--conda-prefix /path/to/build_conda_env`` on snakemake command line.


HPC specifications
==================

To install CulebrONT on a global way, we recommend to charge required dependencies from CulebrONT and modifying *cluster_config.yaml* file.

Preparing *cluster_config.yaml*
-------------------------------

On ``cluster_config.yaml`` , you can add partition, memory and threads to be used by default for each rule. If more memory or threads are requested, please adapt the content of this file before running on a cluster for every rule.

.. warning::
    please adapt the content of this file before running on a cluster for every rule !!

Here is a example of the configuration file we used on the i-Trop HPC.

.. code-block:: yaml

   __default__:
       cpus-per-task : 4
       ntasks : 1
       mem-per-cpu : '2'
       partition : "normal"
       output : 'logs/stdout/{rule}/{wildcards}'
       error : 'logs/error/{rule}/{wildcards}'

   run_nanopolish :
       cpus-per-task : 12
       mem-per-cpu : '4'
       partition : "long"

   run_canu:
       cpus-per-task : 8
       mem-per-cpu : '8'
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
