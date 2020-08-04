Installation
------------

Mandatory installation
----------------------

CulebrONT uses |PythonVersions| |SnakemakeVersions| |Singularity| could be of help (see `below <#singularity>`_\ )


* Install or update CulebrONT

.. code-block:: bash

   git clone https://github.com/SouthGreenPlatform/CulebrONT_pipeline.git
   cd CulebrONT_pipeline

Optional Installation
---------------------

To obtain a clear and correct report, please add also the following dependencies from R:


* ``remote::install_github("strengejacke/strengejacke")``
* 'plotly', 'dplyr', 'optparse', 'htmltools', 'rmdformats', 'magrittr', 'yaml', 'png', 'here', 'htmlwidgets'.

Available data test
-------------------

A data test ``Data-Xoo-sub/`` is available on https://itrop.ird.fr/culebront_utilities/. Feel free to download it using ``wget`` and put it on CulebrONT repertory.


.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7%2B-blue
   :target: https://www.python.org/downloads
   :alt: Python 3.7+

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake 5.10.0+

.. |Singularity| image:: https://img.shields.io/badge/singularity-available-7E4C74.svg
   :target: https://sylabs.io/docs/
   :alt: Singularity 3.10.0+
