Parameters
====================================

A full list of parameters can be found in the table at the bottom of this page. However, in practice, only a few parameters will be relevant for most users of :code:`EUKulele`. These are the required ones:

- mets_or_mags: Whether the user intends to run the analysis for metatranscriptomic samples ("mets") or metagenomic samples ("mags")


.. list-table:: Full list of ``EUKulele`` parameters
   :widths: 25 25 50
   :header-rows: 1

   * - Flag
     - Configuration File Entry
     - Meaning
   * - ``--config``
     - N/A 
     - The path to the configuration file which should be used to retrieve the equivalent of command-line arguments.
   * - ``-m/--mets_or_mags`` 
     - ``mets_or_mags`` 
     - A required flag to indicate whether metatranscriptomic ("mets") or metagenomic ("mags") samples are being used as input.
   * - ``-s/--sample_dir`` 
     - samples 
     - A required flag to indicate where the samples (metagenomic or metatranscriptomic, depending on "mets_or_mags" flag) are located. 
   * - ``-o/--out_dir`` 
     - output 
     - The path to the directory where output will be stored. Defaults to a folder called ``output`` in the present working directory.
   * - ``--reference_dir`` 
     - reference 
     - A flag to indicate where the reference FASTA is stored, or a keyword argument for the dataset to be downloaded and used. Only used if not downloading automatically.
   * - ``--ref_fasta`` 
     - ref_fasta 
     - The name of the reference FASTA file in ``reference_dir``; defaults to reference.pep.fa if not specified, or is set according to the downloaded file if using a keyword argument.
   * - ``--database`` 
     - database 
     - An optional additional argument for specifying the database name. If the database specified is one of the supported databases (currently, "mmetsp", "eukprot", or "phylodb", it will be downloaded automatically. Otherwise, MMETSP is used as a default. 
   * - ``--run_transdecoder``
     - run_transdecoder (set to 0 or 1)
     - An argument for the user to specify whether or not TransDecoder should be used to translate input nucleotide sequences, prior to ``blastp`` being used (i.e., the equivalent protein-protein alignment with the tool of choice). If included in command line or set to 1 in configuration file, ``TransDecoder`` is run. Otherwise, ``blastp`` is run if protein files are found (according to files in the sample directory ending in ``--p_ext`` (below), or ``blastx`` is run if only nucleotide format files are found. 
   * - ``--nucleotide_extension/--n_ext`` 
     - nucleotide_extension 
     - The file extension for samples in nucleotide format (metatranscriptomes). Defaults to .fasta.
   * - ``--protein_extension/--p_ext`` 
     - protein_extension 
     - The file extension for samples in protein format (metatranscriptomes). Defaults to .faa.
   * - ``-f/--force_rerun`` 
     - force_rerun 
     - If included in a command line argument or set to 1 in a configuration file, this argument forces all steps to be re-run, regardless of whether output is already present.
   * - ``--use_salmon_counts`` 
     - use_salmon_counts 
     - If included in a command line argument or set to 1 in a configuration file, this argument causes classifications to be made based both on number of classified transcripts and by counts.
   * - ``--salmon_dir`` 
     - salmon_dir 
     - If ``--use_salmon_counts`` is true, this must be specified, which is the directory location of the ``salmon`` output/quantification files.
   * - ``--names_to_reads`` 
     - names_to_reads 
     - A file that creates a correspondence between each transcript name and the number of ``salmon``-quantified reads. Can be generated manually via the ``names_to_reads.py`` script, or will be generated automatically if it does not exist. \
   * - ``--transdecoder_orfsize`` 
     - transdecoder_orfsize 
     - The minimum cutoff size for an open reading frame (ORF) detected by ``TransDecoder``. Only relevant if ``--use_transdecoder`` is specified.
   * - ``--alignment_choice`` 
     - alignment_choice 
     - A choice of aligner to use, currently ``BLAST`` or ``DIAMOND``.
   * - ``--cutoff_file`` 
     - cutoff_file 
     - A ``YAML`` file, provided in ``src/EUKulele/static/``, that contains the percent identity cutoffs for various taxonomic classifications. Any path may be provided here to a user-specified file.
   * - ``--filter_metric`` 
     - filter_metric 
     - Either evalue, pid, or bitscore (default evalue) - the metric to be used to filter hits based on their quality prior to taxonomic estimation. 
   * - ``--consensus_cutoff`` 
     - consensus_cutoff 
     - The value to be used to decide whether enough of the taxonomic matches are identical to overlook a discrepancy in classification based on hits associated with a contig. Defaults to 0.75 (75%). 
   * - ``--busco_file`` 
     - busco_file 
     - Overrides specific organism and taxonomy parameters (next two entries below) in favor of a tab-separated file containing each organism/group of interest and the taxonomic level of the query. \
   * - ``--organisms``
     - organisms
     - A list of organisms/groups to test the BUSCO completeness of matching contigs for.
   * - ``--taxonomy_organisms`` 
     - taxonomy_organisms 
     - The taxonomic level of the groupings indicated in the list of ``--organisms``; also a list.
   * - ``--individual_or_summary / -i``
     - individual_or_summary 
     - Defaults to summary. Whether BUSCO assessment should just be performed for the top organism matches, or whether the list of organisms + their taxonomies or BUSCO file (above parameters) should be used (individual). When ``-i`` is specified, individual mode is chosen.
   * - ``--busco_threshold``
     - busco_threshold 
     - The threshold for BUSCO completeness for a set of contigs to be considered reasonably BUSCO-complete.
   * - ``--tax_table`` 
     - tax_table 
     - The name of the formatted taxonomy table; defaults to "tax-table.txt.". If this file is not found, it can be generated from the reference FASTA and original taxonomy file using the provided script ``create_protein_file.py``, or the database specified will be automatically downloaded, if it is one of the supported databases.
   * - ``--protein_map``
     - protein_map 
     - The name of the JSON file containing protein correspondences; defaults to "protein-map.json". If this file is not found, it can be generated from the reference FASTA and original taxonomy file using the provided script ``create_protein_file.py``, or the database specified will be automatically downloaded, if it is one of the supported databases.
     