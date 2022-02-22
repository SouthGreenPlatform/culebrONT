.. image:: https://raw.githubusercontent.com/SouthGreenPlatform/culebrONT/master/culebrONT/culebront_logo.png
   :alt: Culebront Logo
   :align: center


|PythonVersions| |SnakemakeVersions| |Singularity|  |Downloads|

.. contents:: Table of Contents
    :depth: 2

About CulebrONT
===============

|readthedocs|

**Homepage:** `https://culebront-pipeline.readthedocs.io/en/latest/ <https://culebront-pipeline.readthedocs.io/en/latest/>`_

Using data from long reads obtained by Oxford Nanopore Technologies
sequencing makes genome assembly easier, in particular to solve repeats
and structural variants, in prokaryotic as well as in eukaryotic
genomes, resulting in increased contiguity and accuracy.

Bunch of softwares and tools are released or updated every week, and a
lot of species see their genome assembled using those.

That’s right.

“*But which assembly tool could give the best results for my favorite
organism?*”

**CulebrONT can help you!** CulebrONT is an open-source, scalable,
modulable and traceable snakemake pipeline, able to launch multiple
assembly tools in parallel and providing help for choosing the best
possible assembly between all possibilities.

Citation
________

https://doi.org/10.1101/2021.07.19.452922

Authors
_______

* Julie Orjuela (IRD)
* Aurore Comte(IRD)
* Sébastien Ravel(CIRAD),
* Florian Charriat(INRAE)
* Sébastien Cunnac(IRD).
* Tram Vi (IRD, AGI)
* Francois Sabot(IRD)

Useful notes
============

Before launching CulebrONT, you could base-calling of arbitrarily
multiplexed libraries across several Minion runs with sequencing quality
control and gather the output files by genome for subsequent steps. For
that use https://github.com/vibaotram/baseDmux.

Thanks
======

Thanks to Ndomassi Tando (i-Trop IRD) by administration support.

The authors acknowledge the `IRD i-Trop HPC <https://bioinfo.ird.fr/>`_ (`South Green Platform <http://www.southgreen.fr>`_) at IRD
Montpellier for providing HPC resources that have contributed to this work.

Thanks to `Yann Delorme <https://nimarell.github.io/resume>`_ for this beautiful logo

License
=======

Licencied under `CeCill-C <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html>`_ and GPLv3.

Intellectual property belongs to `IRD <https://www.ird.fr>`_ , `CIRAD <https://www.cirad.fr/>`_ and authors.

.. |PythonVersions| image:: https://img.shields.io/badge/python-≥3.6%2B-blue
   :target: https://www.python.org/downloads

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥6.10.0-brightgreen.svg
   :target: https://snakemake.readthedocs.io

.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/

.. |readthedocs| image:: https://pbs.twimg.com/media/E5oBxcRXoAEBSp1.png
   :target: https://culebront-pipeline.readthedocs.io/en/latest/
   :width: 400px

.. |Downloads| image:: https://img.shields.io/pypi/dm/culebrONT?color=purple&logo=culebrONT-pypi
   :target: https://pypi.org/project/culebrONT
   :alt: PyPi downloads
