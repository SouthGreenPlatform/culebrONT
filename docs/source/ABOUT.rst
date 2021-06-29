.. image:: _images/culebront_logo.png
   :target: _images/culebront_logo.png
   :align: center
   :alt: Culebront Logo

Today, assembling a genome using long reads from Oxford Nanopore Technologies is really interesting, in particular to solve repeats and structural variants in prokaryotic as well as in eukaryotic genomes. Due to the read length and quality increases, assemblies are increasing themeselves in term of contiguity and accuracy.

We are confronted to a daily deluge of data sequences and of new tools being released or updated every week, allowing many species and samples to have their genomes assembled and that is great ...

"\ *But which assembly tool could give the best results for your favorite organism?*\ "

**CulebrONT can help you!** CulebrONT is an open-source, scalable, modular and traceable Snakemake pipeline, able to launch multiple assembly tools in parallel, giving you the possibility of circularise, polish, correct, and checking quality for assemblies. CulebrONT can help choosing the best assembly pipeline between all possibilities.
Assembly, circularisation, polishing and correction steps are included in CulebrONT, and can be activated (or not) according to user’s requests. The most commonly used tools in the community for each step are integrated, as well as various quality control tools. CulebrONT also generates a report compiling information obtained at each step, to help user to decide which assembly pipeline results to work with.


.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

From assembly to correction
---------------------------

CulebrONT is a really flexible tool to assemble, circularise (or not), polish and correct assemblies. You can give parameters on the *config.yaml* file to CulebrONT to generate a modular, dedicated pipeline for your own data.

.. warning::
    Logically, you must launch at least one of assemblers included in CulebrONT and pipe assembly with circularisation, polishing and correction steps as well as with the quality control pipeline.

Assembly
........

CulebrONT includes (at the moment) six recent and community-validated assemblers.

Included tools :

* Flye version >= 2.6 https://github.com/fenderglass/Flye
* Canu version >= 2.0 https://canu.readthedocs.io/en/latest/quick-start.html
* Miniasm version >= 0.3 https://github.com/lh3/miniasm + Minipolish version >= 0.1.2 https://github.com/rrwick/Minipolish
* Shasta version 0.5.1 https://github.com/chanzuckerberg/shasta
* Smartdenovo version 1.0.0 https://github.com/ruanjue/smartdenovo
* Raven version >= 1.2.2 https://github.com/lbcb-sci/raven


Optional circularisation
........................

Using CulebrONT you can activate or deactivate circularisation steps. Typically, if you are interested on eukaryotic organims, circularisation is not necessary, use CIRCULAR=False on ``config.yaml``  file. Option circular on the configuration file is key on the workflow framework.

Directed acyclic graphs (DAGs) show the differences between deactivated (CIRCULAR=False):

.. image:: _images/schema_pipeline_global-NOCIRC.png
   :target: _images/schema_pipeline_global-NOCIRC.png
   :alt: CIRCULAR_FALSE

and activated CIRCULAR step on configuration file (CIRCULAR=True):

.. image:: _images/schema_pipeline_global-CIRC.png
   :target: _images/schema_pipeline_global-CIRC.png
   :alt: CIRCULAR_TRUE


If an assembled molecule is circular, e.g. for bacteria, this molecule is tagged and will be treated specially in pipeline. We implemented tagging and rotation of circular molecule before each polishing step, and we fixing start position on circular genome. This is efficient when multiple genome alignments are envisaged.

.. note::

    * If circular step is activated, the *--plasmids* option on Flye is used.
    * *Circlator* tool is used to circularise assemblies from Canu and Smartdenovo. Circlator will attempt to identify each circular sequence and output a linearised version from each of them.
    * On circularisation step performed by circlator, trimmed corrected fastq files obtained by Canu are used by circlator. Raw fastq files are used directly for others assemblers.
    * Circularisation for Miniasm is performed by minipolish.
    * Circularisation for Miniasm, Raven and Shasta is checked using generated GFA files on a special tag_circular step, tagging circular fasta sequences.


**OPTIONAL FIXING START**

Only if the circular step is activated, a fixstart step is performed before the quality control. This step uses circlator tool with the option fixstart to rotate circular sequences. Rotation is done using the start dnaA gene (if found). This is important when multiple alignments are envisaged.

Fixstart is always performed on the last draft sequenced obtained on the pipeline. In others words:

.. note::
    * Fixstart is launched after the assembly step if only assembly is activated.
    * Fixstart is launched after the polishing step if only assembly+polishing are activated.
    * Fixstart is launched after correction if assembly+polishing+correction or assembly+correction are activated

.. warning::
    In any case, Fixstart will be deactivated if CIRCULAR is False


Included tools :

* Circlator version >= 1.5.2 https://github.com/sanger-pathogens/circlator


Polishing
.........

Polishing step is ensured by Racon. Racon corrects raw contigs generated by rapid assembly methods with original ONT reads. Choose how many rounds of Racon you want to perform (constrains from 1 to 9 rounds), and CulebrONT will recursively do it for you. Generally 2 to 4 iterations are the best choices.

.. note::

    * Minipolish includes racon polishing, so if polishing is deactivated for the others assemblers (flye, shasta, raven, smartdenovo and canu) miniasm will polish anyway, please take it into account to comparisons.

    * Raven parameter -p (for polishing) is by default 0, racon is performed on CulebrONT to control rotation of circular molecule before every racon step.

Included tools :

* Racon version >= 1.4.13 https://github.com/isovic/racon


Correction
..........

Correction can improve the consensus sequence for a draft genome assembly. We include *Nanopolish* and *Medaka* on correction steps. With CulebrONT you can train a Medaka model and use it directly to obtain a consensus from you favorite organism.

.. note::
    * We have included a split on segments of the assembled molecule before nanopolish and medaka. Each segment is polished on parallel to improve speed and gain time. Segments polished are merged subsequently. CulebrONT has implemented parallelism following `medaka documentation <https://nanoporetech.github.io/medaka/installation.html#improving-parallelism>`_ and `nanopolish practices <https://nanopolish.readthedocs.io/en/latest/quickstart_consensus.html#compute-a-new-consensus-sequence-for-a-draft-assembly>`_.

If you have short reads, you can now use *Pilon* to correct assemblies. As racon, using CulebrONT several recursive rounds of Pilon can now be run !

Included tools :

* Medaka Medaka-gpu version >= 1.2 https://github.com/nanoporetech/medaka
* Nanopolish version >= 0.13.2 https://nanopolish.readthedocs.io/en/latest/index.html#
* Pilon version >= 1.24 https://github.com/broadinstitute/pilon/releases/


Quality on assemblies
----------------------

A variety of useful tools exist to check high accuracy assemblies.

.. image:: _images/schema_pipeline_global-QUALITY.png
   :target: _images/schema_pipeline_global-QUALITY.png
   :alt: QUALITY


CulebrONT checks the quality of the assemblies with using these optional tools:

.. note::

    * BUSCO: helps to check if you have a good assembly, by searching the expected single-copy lineage-conserved orthologs in any newly-sequenced genome from an appropriate phylogenetic clade.
    * QUAST: is a good starting point to evaluate the quality of assemblies, providing many helpful contiguity statistics.
    * Blobtools: allows to detect contamination on assembled contigs as well as CG% biases.
    * Assemblytics: compares structural variations of assemblies versus a reference genome
    * KAT: explores kmers frequencies and checks possible contamination
    * Samtools flagstats: calculates remapping stats using illumina reads over assemblies
    * Mauve: allows multiple alignment of several assembles (for small genome only).


.. danger::

    Please, activate Mauve only for small genomes.

Included tools :

* BUSCO version >= 4.0.5
* QUAST version >= 5.0.2
* Bloobtools version >= 1.1.1
* Assemblytics version >= 1.2
* KAT version >= 2.4.2
* Samtools version>= 1.10
* Mauve > 2.4.0.snapshot_2015_02_13
