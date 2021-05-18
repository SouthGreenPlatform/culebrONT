
![Culebront Logo](./docs/source/SupplementaryFiles/culebront_logo.png)


[![PythonVersions](https://img.shields.io/badge/python-3.7%2B-blue)](https://www.python.org/downloads)
[![SnakemakeVersions](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Singularity](https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg)](https://sylabs.io/docs/)
[![Conda](https://img.shields.io/badge/conda-4.8.5%20-green)](https://docs.conda.io/projects/conda/en/latest/index.html)

Using data from long reads obtained by Oxford Nanopore Technologies sequencing makes genome assembly easier, in particular to solve repeats and structural variants, in prokaryotic as well as in eukaryotic genomes, resulting in increased contiguity and accuracy.

Bunch of softwares and tools are released or updated every week, and a lot of species see their genome assembled using those.

That’s right.

"*But which assembly tool could give the best results for my favorite organism?*"

**CulebrONT can help you!** CulebrONT is an open-source, scalable, modulable and traceable snakemake pipeline, able to launch multiple assembly tools in parallel and providing help for choosing the best possible assembly between all possibilities.

**Homepage: https://culebront-pipeline.readthedocs.io/en/latest/**

<a name="citation"></a>
## Citation

@Authors:

Julie Orjuela (IRD), Aurore Comte(IRD), Sébastien Ravel(CIRAD), Florian Charriat(INRAE), Tram Vi(IRD, AGI), Francois Sabot(IRD) and Sébastien Cunnac(IRD).

<a name="notes"></a>

## Useful notes

Before launching CulebrONT, you could base-calling of arbitrarily multiplexed libraries across several Minion runs with sequencing quality control and gather the output files by genome for subsequent steps. For that use https://github.com/vibaotram/baseDmux.

#### Thanks

Thanks to Ndomassi Tando (i-Trop IRD) by administration support.

The authors acknowledge the IRD i-Trop HPC (South Green Platform) at IRD Montpellier for providing HPC resources that have contributed to this work. https://bioinfo.ird.fr/ - http://www.southgreen.fr

Thanks to Yann Delorme for this beautiful logo https://nimarell.github.io/resume

<a name="licence"></a>
## License
Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
Intellectual property belongs to IRD and authors.
