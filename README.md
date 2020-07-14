# EUKulele

![EUKulele-logo](/.infrastructure/eukulele-logo.png)

## Formalizing environmental eukaryotic taxonomic assignment

### About EUKulele
`EUKulele` is a Python program for taxonomic annotation of microbes in metatranscriptomic and metagenomic samples, with special emphasis on eukaryote discovery. `EUKulele` can be downloaded from [PyPI](https://pypi.org/) and used as a Python module, or it may be downloaded via `conda` and used as a command-line program. The software includes three major features:
- Database setup and formatting
- Database creation, alignment, and taxonomic estimation
- Assessment of the BUSCO completeness of subsets of contigs at each taxonomic level

### Prerequisites for running EUKulele
In principle, there are two prerequisites for running the software:
1. Metagenomic or metatranscriptomic sample files (unless using the provided sample data)
2. A database to align the contigs from the metagenome/metatranscriptome to

Three databases are supported by default from within `EUKulele`, and may be downloaded and formatted automatically if the user chooses:
- [PhyloDB](https://drive.google.com/drive/u/0/folders/0B-BsLZUMHrDQfldGeDRIUHNZMEREY0g3ekpEZFhrTDlQSjQtbm5heC1QX2V6TUxBeFlOejQ)
- [EukProt](https://figshare.com/articles/EukProt_a_database_of_genome-scale_predicted_proteins_across_the_diversity_of_eukaryotic_life/12417881/2)
- [MMETSP](https://zenodo.org/record/1212585#.Xw3PoJNKhTZ)

### EUKulele parameters
Parameters may either be passed as arguments directly to EUKulele, or may be specified in a configuration `YAML` file. If running `EUKulele` as a Python package, you can call the function `eukulele` with two potential arguments, one of which is required:

'''
eukulele(config="", string_arguments="")
'''

Where the `config` argument allows you to list the path of the configuration `YAML` file that you wish to use, or the `string_arguments` option allows you to create a string containing all of the arguments that you would use normally using the command-line options specified below. 

`EUKulele` can be run from the command line using: `EUKulele [arguments]`, once installed. 

## Previous README file notes

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
