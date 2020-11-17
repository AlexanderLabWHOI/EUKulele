Expected Output of ``EUKulele``
================================

Below is what you should expect to see when you run ``EUKulele``. ``output-folder-name`` is the top level directory specified by the ``-o/--out_dir`` parameter, which defaults to a folder called "output" in the current working directory.

| output-folder-name
| ├── busco
| │   ├── build
| │   ├── make.bat
| │   ├── Makefile
| │   └── source
| ├── busco_assessment
| ├── core_taxonomy_counts
| ├── core_taxonomy_estimation
| ├── core_taxonomy_visualization
| ├── mets (if running ``TransDecoder``, or protein files provided)
| │   ├── *.faa files*
| │   └── transdecoder
| │   |   └── *other TransDecoder outputs*
| ├── mets_full / mags_full (depending on ``--mets_or_mags`` parameter)
| │   └── diamond (or blast)
| ├── mets_core / mags_core (depending on ``--mets_or_mags`` parameter)
| │   └── diamond (or blast)
| ├── taxonomy_counts
| ├── taxonomy_estimation
| └── taxonomy_visualization
|


Taxonomy Estimation Folders
---------------------------

Inside each of the taxonomy estimation folders (``core_taxonomy_estimation``, for exclusively transcripts annotated as core genes, and ``taxonomy_estimation``), there are files labeled ``<sample_name>-estimated-taxonomy.out``. Each of these files has the following columns:

- transcript_name
    - The name of the matched transcript/contig from this sample file
- classification_level
    - The most specific taxonomic level that the transcript/contig was classified at
- full_classification
    - The full taxonomic classification for the transcript/contig, up to and including the most specific classification
- classification
    - Just the most specific classification obtained, at the taxonomic level specified by classification_level
- max_pid
    - The maximum percentage identity reported by the alignment program for the match
- counts
    - 1 if not using ``Salmon``, or the number of counts reported by the quantification program for that transcript/contig if using ``Salmon``
- ambiguous
    - 1 if there were multiple disagreeing matches for this transcript/contig, resolved by either consensus annotation using the user-provided cutoff (defaults to 75%) or by last common ancestor (LCA) approaches, otherwise 0

Taxonomy Counts Folders
-----------------------

Inside each of the taxonomy counts folders (``core_taxonomy_counts`` and ``taxonomy_counts``), there is a file with the aggregated counts belonging to each annotation at the possible taxonomic levels. These files are derived diectly from the taxonomy estimation and used for the visualization steps. Inside this folder are comma-separated files with the following headings:

- *Taxonomic Level*
    - The taxonomic level specified by the file name
- GroupedTranscripts
    - A semicolon-separated list of the identification names for transcripts that matched to this taxonomic label
- NumTranscripts
    - The total number of transcripts (the length of the semicolon-separated list of transcripts)
- Counts
    - Provided if ``Salmon`` counts are used - the matching summed number of counts for each label
- Sample
    - The original metagenomic/metatranscriptomic sample that this count is from (a separate row would be provided if the match is found in multiple samples)
    
The taxonomic count files are named according to the convention ``<output-folder-name>_all_<taxonomic-level>_counts.csv``.

Taxonomy Visualization Folders
------------------------------

Inside each of the taxonomy visualization folders (``core_taxonomy_visualization``, for exclusively transcripts annotated as core genes, and ``taxonomy_visualization``), there are auto-generated barplots that show:

- x-axis: samples
- y-axis, left subplot (if using counts): relative number of transcripts
- bars, left subplot (if using counts): each of the top represented taxonomic groups (must represent >= 5% of total number of transcripts)
- y-axis, right subplot (if using counts): relative number of counts
- bars, right subplot (if using counts): each of the top represented taxonomic groups (must represent >= 5% of total counts)

The right subplot is only generated if counts from a quantification tool (namely, ``Salmon``) are provided.
