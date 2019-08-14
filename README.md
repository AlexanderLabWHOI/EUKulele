# EUKulele
## Formalizing environmental eukaryotic taxonomic assignment

## General flow

- Users input .fa files into `samples/` directory. Samples can be either metaT or metaG bins. Users select which type in the config files.
- Users also specify a database that they want to annotate against. Database needs to be setup such that there are unique identifiers for each transcriptome, genome, etc. that is included. Files should be provided in... (nucleotide? protein?? Leaning towards protein.)
- Associated with the reference database there should be taxon flow with the format eukarota; xxx; xxx;  in a standardized format.
  - We will provide a zenodo download link that has a pre-formatted database of eukaryotes with the specified file formats.
  - We could/should create a helper function / script to automate the generation of this file type?
  - We might consider looking into BIOM files but this might be overkill at this point. Helps integrate with an ecosystem of other tools though (e.g. qiime2)
- For MAGs:
  - For every predicted protein calculate LCA based on taxonomy % hit (% identity to DB hit)
  - Amalgamate all the hits for the proteins and assess the LCA based on all hits. % cutoffs? How to do this?
  - Output all proteins with top hits and proteins that do not have good matches.
  - Provide stats on next best hit etc.
- For metaT:
  - For every contig (transcript) predict protein --> run diamond alingment against reference DB
  - Estimate LCA based on % id
  - Group cotigs/transcripts into variuos functional units. Allow users to specify groupings?
