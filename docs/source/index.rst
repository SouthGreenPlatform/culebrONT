Welcome to CulebrONT's documentation!
=====================================


.. image:: _images/culebront_logo.png
   :target: _images/culebront_logo.png
   :align: center
   :alt: Culebront Logo

Today, assembly a genome using long reads from either Oxford Nanopore Technologies or Pacific Biosciences is really powerful, in particular to solve repeats and structural variants, for prokaryotic as well as for eukaryotic genomes. Such technologies provide assemblies that are increased in contiguity and accuracy.

Due to the daily deluge of data sequences and the increasing number of released tools that are even updated every week, many species see having their genome assembled in almost chromosome-scale, and that’s great...

However a huge question remains:

"\ *But which assembly tool will provide the best result for your favorite organism?*\ "

To that anguishing idea, we can answer: **CulebrONT can help you!**

CulebrONT is an open-source, scalable, modular and traceable *Snakemake* pipeline, able to launch multiple assembly tools in parallel, giving you the possibility of circularise, polish, and correct assemblies, in addition to perform quality controls. CulebrONT can help to choose the best assembly pipeline between all possibilities.

.. toctree::
   :caption: About CulebrONT
   :name: about_culebront
   :maxdepth: 2

   ABOUT.rst

.. toctree::
   :caption: Install
   :name: install
   :maxdepth: 2

   INSTALL.rst

.. toctree::
   :caption: Defining Workflows
   :name: defining_workflows
   :maxdepth: 2

   WORKFLOWS.rst

.. toctree::
   :caption: Project Info
   :name: project_info
   :maxdepth: 1

   PROJECT.rst

.. toctree::
   :caption: FAQ
   :name: faq_questions
   :maxdepth: 2

   FAQ.rst

.. toctree::
   :caption: API
   :name: api_culebront
   :maxdepth: 2

   click.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
