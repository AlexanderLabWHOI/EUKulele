.. _usingeukulele:
====================================
Using EUKulele
====================================

``EUKulele`` has been designed to provide taxonomic annotation for two primary data types: 1) contigs derived from metatranscriptomes (METs) and 2) metagenome assembled genomes (MAGs). The focus of ``EUKulele`` is on the annotation of microbial eukaryotes; however, it is conceivable to use any database as the foundation of your analyses (see Databases).

Metatranscriptomes (METs)
=========================
In the first case, metatranscriptomes (shortened in ``EUKulele`` to ``mets``), are assumed to be contigs generated from shotgun-style sequencing and assembly of metatranscriptomic data (RNA) from a mixed community. These contigs can be provide to ``EUKulele`` as either nucleotide sequences (such as those output by `Trinity <https://github.com/trinityrnaseq/trinityrnaseq/wiki>`_) or predicted protein sequences from these contigs (such as those output by `Transdecoder <https://github.com/transdecoder>`_). 

The most basic running of ``EUKulele`` on metatranscriptome samples  would be::

    EUKulele --sample_dir path/to/metatranscriptome/samples -m mets

This command will automatically download (if not downloaded already) the MMETSP database and use ``DIAMOND`` to align sequences to the database. ``EUKulele`` will automatically inspect the directory for files ending with ``.fna`` (nucleotide files) or ``.faa`` (amino acid files). If no files are found within the directory with those extensions an error will be raised. If nucleotide files are found (``.fna``) a blastx-style search will automatically performed; alternatively, if amino acid files are found a blastp-style search will be automatically performed.  The results of the ``DIAMOND`` searches and estimated taxonomy for each of the contigs within the input datasets as well as various figures will be output in a directory called ``output``. More on output files in the :ref:`Output Documentation Section<documentation>`! 

``EUKulele`` is highly customizable and can be easily adapted to work with your dataset. For example, if you wanted to run ``BLAST`` against the ``PhyloDB`` database on your protein files that end in ``.protein`` this can be accomplished with the following command::

    EUKulele --sample_dir path/to/metatranscriptome/samples -m mets --protein_extension .protein --database phylodb --alignment_choice BLAST

A full list of parameters and customizations can be found here :ref:`Parameters Section<parameters>`.  Users might find it simpler  to modify a config file such as that provided `here <https://github.com/AlexanderLabWHOI/EUKulele/blob/master/config.yaml>`_). A config file can be passed to ``EUKulele`` as follows. A config file will be automatically generated after the first run of ``EUKulele`` in you working directory to facilitate re-running:: 

    EUKulele --config config.yaml

.. note::
    It is feasible to run ``EUKulele`` on metagenome derived contigs (not MAGs, those are discussed in detail below). If you wish to analyze metagenomic contigs as is described above for metatranscriptomes, we **STRONGLY RECOMMEND** that you provide predicted proteins from your metagenome rather than the nucleotide sequences from your metagenomic assembly. 
    Metagenomic contigs often consist of many proteins. ``EUKulele`` can predict proteins from metatranscriptomes with ``TransDecoder``, but this is **NOT** advised for metagenomic contigs. Additionally, a ``blastx``-style search will no be optimal for metagenomic contigs. It is up to the user to provide predicted proteins from their metagenomic contigs, as this can be a complex process (particularly for eukaryotic metagenomes) and is not within the scope of ``EUKulele``.
    Note that this is a problem in particular for eukaryotic sequences, where protein calling is more complicated due to the presence of introns.

Metagenome Assembled Genomes (MAGs)
===================================
In the second case, metagenome assembled genomes (shortened in ``EUKulele`` to ``mags``), are assumed to be the predicted proteins from a MAG derived from the binning of like contigs from a metagenome assembly based on tetranucleotide frequency, abundance, etc. as done by programs like `CONCOCT <https://github.com/BinPro/CONCOCT>`_ or `metaBAT <https://bitbucket.org/berkeleylab/metabat>`_. Please note that we only recommend using predicted proteins for the taxonomic annotation of MAGs and not raw nucleotide sequences (see above note for some justification). 

The most basic invocation of a MAG based analysis with EUKulele is the following::

    EUKulele --sample_dir path/to/MAGs -m mags

This will align all MAG proteins against the ``MMETSP`` with ``DIAMOND``. As detailed above in the section on metatranscriptomes, the parameters (see :ref:`Parameters Section<parameters>`) or config file can be used to adjust the running of ``EUKulele`` to suit the user. 

Of particular interest to the user might be the ``--tax_cutoffs`` parameter and associated file (``tax-cutoffs.yaml``). The default ``tax-cutoffs.yaml`` file details the percent id cutoffs that are used to assign taxonomy at various levels::

    species: 95
    genus: 80
    family: 65
    order: 50
    class: 30

These cutoffs are based on commonly used cutoffs in the literature, but can be modified by the user to suit their purposes. Using the above example, proteins that align with 95% identity to a reference protein will be assigned at the level of the species, those aligning with 65% identity at the level of family, etc. 

As with the metatranscriptome analysis, the taxonomic estimation of every protein will be output in the ``taxonomy_estimation`` folder within the output directory. Additionally, however, the taxonomic consensus across all proteins within a MAG will also be returned for MAGs. 

.. list-table:: Example ``max-level-mag`` output file
   :widths: 25 25 25
   :header-rows: 1

   * - 
     - max_taxa
     - proportion_id
   * - supergroup 
     - Eukaryota      
     -  0.9998
   * - division 
     - Archaeplastida      
     - 0.9898
   * - class   
     - Chlorophyta  
     - 0.9878
   * - family   
     - Mamiellales  
     - 0.9764
   * - genus   
     - Mamiellaceae  
     - 0.9382
   * - species   
     - Micromonas sp. NEPCC29
     - 0.4002

The ``max-level-mag`` files detail the relative proportion of all the proteins within a MAG that agree at a particular level. The file reports each of the six level considered by default in ``EUKulele`` (more on this in the section on databases :ref:`here<databases>`). For each level, the taxa which recruited that most reads within that level is reported. For example, the majority of proteins in the division level were annotated as Archaeplastida. The proportion of proteins that are annotated as that max level are also reported. 

So, in the above example 99.98% of the proteins in the dataset have a best hit to the supergroup level Eukaryota, meaning that the vast majority of the proteins had the same annotation at the supergroup level. This is largely true, where all proteins are annotated consistently (>90%) from supergroup to genus. However, only 40% of the proteins annotated consistently at the species level. It is up to the user to decide where and how they want to make a final taxonomic annotation for their MAG. In the above example, one might choose to annotate with confidence to the level of genus given the universally high consensus across proteins.

LCA Algorithm
=============

In some cases, multiple hits from alignment via ``blast`` or ``diamond`` will be reported and will meet the threshold specified by the user (see the :ref:`Parameters Section<parameters>`). In this case, the hits available at each taxonomic level will be evaluated using a simple Last Common Ancestor (LCA) algorithm. This simple implementation of the algorithm accepts input from the user (detailed in the :ref:`Parameters Section<parameters>`; parameter is ``--consensus_cutoff`` and has default of 0.75/75%) on what percentage of alignment-derived annotations need to be identical in order for the annotation to be adopted. If, for instance, only 50% of alignment hits match at the species level, less specific taxonomic levels are assessed until a 75% consensus is reached. For example, if two of four hits have the same species annotation, but all four hits have the same genus annotation, the genus annotation would be used, even if all hits meet the defined percentage identity threshold for the species level. 

LCA, while a robust annotation approach, is not the only means of predicting taxonomic level. We are currently exploring adding a phylogenetic estimate of eukaryotic taxonomy, particularly for the taxonomic placement of MAGs.
