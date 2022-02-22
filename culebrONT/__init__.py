#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .snakeWrapper import *
from culebrONT.global_variables import *

logo = INSTALL_PATH.joinpath('culebront_logo.png').as_posix()


__version__ = get_version()

__doc__ = """Today, assembly a genome using long reads from Oxford Nanopore Technologies is really interesting in 
particular to solve repeats and structural variants in prokaryotic as well as in eukaryotic genomes. Assemblies are 
increasing contiguity and accuracy. The daily increase of data sequences obtained and the fact that more and more 
tools are being released or updated every week, many species are having their genomes assembled and that’s is great … 
“But which assembly tool could give the best results for your favorite organism?” CulebrONT can help you! CulebrONT 
is an open-source, scalable, modular and traceable Snakemake pipeline, able to launch multiple assembly tools in 
parallel, giving you the possibility of circularise, polish, and correct assemblies, checking quality. CulebrONT can 
help to choose the best assembly between all possibilities. """

description_tools = f"""
    Welcome to CulebrONT version: {__version__} ! Created on November 2019
    @author: Julie Orjuela (IRD), Aurore Comte (IRD), Sebastien Ravel (CIRAD), Florian Charriat (INRAE),
    Bao Tram Vi (IRD), François Sabot (IRD) and Sebastien Cunnac (IRD)
    @email: julie.orjuela@ird.fr, aurore.comte@ird.fr

    #     .-+.
    #   `odNNh
    #   +NNd:
    #  .Nh.   ---:`
    #  -Nm`  ````./
    #   oNm/ ```-o-
    #    .yNmo/+/.     .oooooo.               oooo             .o8                  .oooooo.   ooooo      ooo ooooooooooooo
    #    `-+ydNmo.    d8P'  `Y8b              `888            "888                 d8P'  `Y8b  `888b.     `8' 8'   888   `8
    #  `/s+../oNNN:  888          oooo  oooo   888   .ooooo.   888oooo.  oooo d8b 888      888  8 `88b.    8       888
    #  +s/     `hNm  888          `888  `888   888  d88' `88b  d88' `88b `888""8P 888      888  8   `88b.  8       888
    #  os:  ////hNN` 888           888   888   888  888ooo888  888   888  888     888      888  8     `88b.8       888
    #  -so- .//sNNs  `88b    ooo   888   888   888  888    .o  888   888  888     `88b    d88'  8       `888       888
    #   `/oo:.yNd/    `Y8bood8P'   `V88V"V8P' o888o `Y8bod8P'  `Y8bod8P' d888b     `Y8bood8P'  o8o        `8      o888o
    #     -yooo+.
    #   `yNs.`-/oo:
    #   dNo` ....+s+
    #   :shmdddhhy+:

    Please cite our github: {GIT_URL}
    Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
    and GPLv3 Intellectual property belongs to IRD, CIRAD and authors.
    Documentation avail at: {DOCS}
    {get_last_version(url=GIT_URL, current_version=__version__)}"""


MODULE_FILE = f"""#%Module1.0
##

## Required internal variables
set     prefix       {INSTALL_PATH.as_posix().strip()}
set     version      {__version__.strip()}

# check if install directory exist
if {{![file exists $prefix]}} {{
    puts stderr "\t[module-info name] Load Error: $prefix does not exist"
    break
    exit 1
}}

## List conflicting modules here
conflict culebrONT

## List prerequisite modules here
module load graphviz/2.40.1

set		fullname	CulebrONT-{__version__.strip()}
set		externalurl	"\n\t{DOCS.strip()}\n"
set		description	"\n\t{__doc__.strip()}

## Required for "module help ..."
proc ModulesHelp {{}} {{
  global description externalurl
  puts stderr "Description - $description"
  puts stderr "More Docs   - $externalurl"
}}

## Required for "module display ..." and SWDB
module-whatis   "loads the [module-info name] environment"

## Software-specific settings exported to user environment

prepend-path PATH $prefix

"""
