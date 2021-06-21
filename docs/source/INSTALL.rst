.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Requirements
============

CulebrONT has been mostly developed to work on an HPC. Let's see how to install it.

CulebrONT requires |PythonVersions|, |SnakemakeVersions| and |graphviz|.

Additionally, as Culebront uses many tools, you must install and use them through two possibilities:

1. Either through the |Singularity| containers,

2. Or using the ``module load`` mode,

3. Or you can also mix both! (yes we know it is a third possibility :D)

Steps for installation
======================

First of all clone our repository or download last version from https://github.com/SouthGreenPlatform/CulebrONT_pipeline.

.. code-block:: bash

   git clone https://github.com/SouthGreenPlatform/CulebrONT_pipeline.git
   cd CulebrONT_pipeline


When you have cloned it and went into the folder,

1. Adapt the file :file:`tools_path.yaml` - in YAML (Yet Another Markup Language) - format  to indicate to CulebrONT where the tools are installed.
See the section :ref:`1. How to configure tools_path.yaml` for details.


2. Adapt the :file:`cluster_config.yaml` file to manage cluster resources such as partition, memory and threads available for each job.
See the section :ref:`2. Preparing *cluster_config.yaml*` for further details.


3. Create a *snakemake profile* to configure cluster options.
See the section :ref:`3. Snakemake profiles` for details.

4. Create CulebrONT envmodules file.
See the section :ref:`4. Export CulebrONT to $PATH` for details.

5. Modify run script. See the section :ref:`5. Adapt `submit_culebront.sh`` for details.

6. Test CulebrONT install (Optional but recommended) using an available dataset.
See the section :ref:`6. Check install` for details.


1. How to configure tools_path.yaml
-----------------------------------

In the ``tools_path.yaml`` file, you can find two sections: SINGULARITY and ENVMODULES. In order to fill it correctly, you have 3 options:

1. Use only SINGULARITY containers. In this case, fill only this section. Put the path to the builded singularity images.
Absolute paths are strongly recommended but not mandatory. See the section :ref:`'How to build singularity images'<How to build singularity images>`  for further details.

.. literalinclude:: ../../tools_path.yaml
    :language: YAML
    :lines: 2-5

.. warning::
    The use of SINGULARITY is constraint with the use the *--use-singularity* parameter in the snakemake command line.

2. Use only ENVMODULES. In this case, fill this section with modules available on your cluster (here is an example):

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


CulebrONT needs a wide set of R modules for reporting, if you use ENVMODULE R. Just have a look at dependencies in the ``Containers/Singularity.report.def`` file.
Yes, plenty of packages!! That's why we provide build singularity containers ready to use and recommend them for the R part.

.. warning::
    The use of ENVMODULE is constraint with the use the *--use-envmodules* parameter in the snakemake command line.
    More details can be found here: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules


3. Use both SINGULARITY and ENVMODULE. In this case, Snakemake will try to use the modules as a priority, and only in a second time singularity containers.
So, you can mix them! It's up to you, depending which modules are available on your favorite cluster.


.. NOTE::

    If you do not need a particular tool, please leave the path empty as shown below for WEESAM.

    .. code-block:: YAML

        SINGULARITY:
            REPORT : './Containers/Singularity.report.sif'
            WEESAM : ''
            TOOLS: './Containers/Singularity.conda.sif'


How to build singularity images
...............................

You can obtain singularity builded images in two ways. Choose one of them (but option 1 is better).

1. Download already builded available singularity images from our i-Trop server https://itrop.ird.fr/culebront_utilities/
It can be long if your internet connexion is bad. Do not forget it stakes hard drive space!

.. code-block:: bash

    cd CulebrONT_pipeline/Containers
    wget --no-check-certificate -rm -nH --cut-dirs=2 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/singularity_build/


2. You can build your own image using the available *.def* recipes from the ``CulebrONT_pipeline/Containers`` directory.

.. warning::
    Be careful, you need root rights to build singularity images

.. code-block:: bash

    cd CulebrONT_pipeline/Containers/
    sudo make build


2. Preparing *cluster_config.yaml*
----------------------------------

In the ``cluster_config.yaml`` file, you can add partition, memory and threads to be used by default or specifically for each rule/tool.


.. warning::
    If more memory or threads are requested, please adapt the content of this file before running on your cluster.


Here is a example of the configuration file we used on our :download:`i-Trop HPC<../../cluster_config.yaml>`.





3. Snakemake profiles
---------------------

The Snakemake-profiles project is an open effort to create configuration profiles allowing to execute Snakemake in various computing environments
(job scheduling systems as Slurm, SGE, Grid middleware, or cloud computing), and available at https://github.com/Snakemake-Profiles/doc.

In order to run CulebrONT on HPC cluster, we strongly recommended to use profiles.

Quickly, here is an example of the Snakemake SLURM profile we use for the French national bioinformatics infrastructure at IFB.
We followed the documentation found here https://github.com/Snakemake-Profiles/slurm#quickstart.


.. code-block:: bash

    $ cd CulebrONT_pipeline/profiles
    $ cookiecutter https://github.com/Snakemake-Profiles/slurm.git
    profile_name [slurm]: CulebrONT
    sbatch_defaults []: --export=ALL --job-name={rule}
    cluster_config []: /path/to/CulebrONT_pipeline/cluster_config.yaml
    Select advanced_argument_conversion:
    1 - no
    2 - yes
    Choose from 1, 2 [1]: 1
    cluster_name []:

Now, your basic profile is created. To finalize it, change the ``CulebrONT_pipeline/profiles/CulebrONT/config.yaml`` to :

.. code-block:: ini

    restart-times: 0
    jobscript: "slurm-jobscript.sh"
    cluster: "slurm-submit.py"
    cluster-status: "slurm-status.py"
    jobs: 200                   # edit to limit number jobs submit
    max-jobs-per-second: 1
    max-status-checks-per-second: 10
    local-cores: 1
    latency-wait: 6000
    use-singularity: True       # if False, please install all R package on tools_config.yaml ENVMODULE/R
    use-envmodules: True        # Adapt True/False
    rerun-incomplete: True
    printshellcmds: True

.. warning::

    If you decided to create a profile in another path than ``CulebrONT_pipeline/profiles/CulebrONT``.
    Don't forget to change also the profile path in ``submit_culebront.sh``.

.. danger::
    If you used profiles with:
    *use-singularity: False* AND *use-envmodules: False*, all CulebrONT tools must be available on $PATH !!!!!!!!!

.. note::

    Note that the *--profile* argument can be either a relative or an absolute path.
    In addition, snakemake will search for a corresponding folder *profile_name* in /etc/xdg/snakemake and in $HOME/.config/snakemake, where globally accessible profiles can be placed.

.. note::
   You can find the generated files  when profiles are created following this documentation. IFB profile example can be found on the `̀ gift_files/CulebrONT_SLURM`` repertory on our github repository as well as a CulebrONT module example to make life easier for our cluster administrators ^^.

4. Export CulebrONT to $PATH
----------------------------

To run script `submit_culebront.sh` you need export path of installation.

.. code-block:: bash

    # used for script submit_culebront.sh
    export CULEBRONT="/path/to/CulebrONT_pipeline"

    # to call submit_culebront.sh script anywhere
    export PATH=$CULEBRONT:$PATH

or you can adapt the module load file :download:`available here<../../gift_files/CulebrONT_envmodules>`

5. Adapt `submit_culebront.sh`
------------------------------

Finally, modify the header of `submit_culebront.sh` file according to your scheduler.

.. literalinclude:: ../../submit_culebront.sh
    :language: sh
    :lines: 1-8

6. Check install
----------------

Optionally, in order to test your install of CulebrONT pipeline, a data test called ``Data-Xoo-sub/`` is available on https://itrop.ird.fr/culebront_utilities/.
Feel free to download it using ``wget`` and put it on CulebrONT directory.

.. code-block:: bash

    cd CulebrONT_pipeline
    wget --no-check-certificat -rm -nH --cut-dirs=1 --reject="index.html*" --no-parent https://itrop.ird.fr/culebront_utilities/Data-Xoo-sub/



Now, it is time to prepare the configuration file ``config.yaml`` file to indicate to CulebrONT which kind of pipeline you want to create and to run on your data !


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
