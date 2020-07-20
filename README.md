# EUKulele

[![MIT License][license-shield]][license-url]
![Python version][python-version-url]

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

### Downloading EUKulele


### EUKulele parameters
Parameters may either be passed as arguments directly to EUKulele, or may be specified in a configuration `YAML` file. If running `EUKulele` as a Python package, you can call the function `eukulele()` with two potential arguments, `config` or `string_arguments`, one of which is required, and then a third argument `command_line`, which should never be set manually and is only relevant if the script has been run from the command line. 

```
eukulele(config="", string_arguments="", )
```

where the `config` argument allows you to list the path of the configuration `YAML` file that you wish to use, while the `string_arguments` option allows you to create a string containing all of the arguments that you would use normally using the command-line options specified below. 

`EUKulele` can be run from the command line using: `EUKulele [arguments]`, once installed. The following arguments are supported, either as command-line arguments, as a string of arguments passed to the `eukulele()` function within the `EUKulele` Python module, or as part of a confirmation `YAML` file, with names specified in the table. A configuration file may also be provided when `EUKulele` is run from the command line, when the `--config` flag is used with a configuration file path.

| Flag 	| Configuration File Entry 	| Meaning 	|
|-	|-	|-	|
| --config 	| N/A 	| The path to the configuration file which should be used to retrieve the equivalent of command-line arguments. 	|
| --mets_or_mags 	| mets_or_mags 	| A required flag to indicate whether metatranscriptomic ("mets") or metagenomic ("mags") samples are being used as input. 	|
| --reference_dir 	| reference 	| A required flag to indicate where the reference FASTA is stored, or a keyword argument for the dataset to be downloaded and used. 	|
| --sample_dir 	| samples 	| A required flag to indicate where the samples (metagenomic or metatranscriptomic, depending on "mets_or_mags" flag) are located. 	|
| --ref_fasta 	| ref_fasta 	| The name of the reference FASTA file in `reference_dir`; defaults to reference.pep.fa if not specified, or is set according to the downloaded file if using a keyword argument. 	|
| --out_dir 	| output 	| The path to the directory where output will be stored. Defaults to a folder called `output` in the present working directory. 	|
| --database 	| database 	| An optional additional argument for specifying the database name. 	|
| --nucleotide_extension/--n_ext 	| nucleotide_extension 	| The file extension for samples in nucleotide format (metatranscriptomes). Defaults to .fasta. 	|
| --protein_extension/--p_ext 	| protein_extension 	| The file extension for samples in protein format (metatranscriptomes). Defaults to .faa. 	|
| --scratch 	| scratch 	| The directory where temporary files will be stored. Defaults to `../scratch`. 	|
| --ref_fasta_ext 	| ref_fasta_ext 	| The file extension for `ref_fasta`, if applicable. Defaults to .fasta. 	|
| --force_rerun 	| force_rerun 	| If included in a command line argument or set to 1 in a configuration file, this argument forces all steps to be re-run, regardless of whether output is already present. 	|
| --use_salmon_counts 	| use_salmon_counts 	| If included in a command line argument or set to 1 in a configuration file, this argument causes classifications to be made based both on number of classified transcripts and by counts. 	|
| --salmon_dir 	| salmon_dir 	| If `--use_salmon_counts` is true, this must be specified, which is the directory location of the `salmon` output/quantification files. 	|
| --names_to_reads 	| names_to_reads 	| A file that creates a correspondence between each transcript name and the number of `salmon`-quantified reads. Can be generated manually via the `names_to_reads.py` script, or will be generated automatically if it does not exist. 	|
| --transdecoder_orfsize 	| transdecoder_orfsize 	| The minimum cutoff size for an open reading frame (ORF) detected by `TransDecoder`. 	|
| -p 	| choose_parallel 	| Can be set to 1 or "parallel" to indicate that taxonomy estimation should be run in parallel. 	|
| --alignment_choice 	| alignment_choice 	| A choice of aligner to use, currently `BLAST` or `DIAMOND`. 	|
| --cutoff_file 	| cutoff_file 	| A `YAML` file, provided in `EUKulele/static/`, that contains the percent identity cutoffs for various taxonomic classifications. Any path may be provided here to a user-specified file. 	|
| --filter_metric 	| filter_metric 	| Either evalue, pid, or bitscore (default evalue) - the metric to be used to filter hits based on their quality prior to taxonomic estimation. 	|
| --consensus_cutoff 	| consensus_cutoff 	| The value to be used to decide whether enough of the taxonomic matches are identical to overlook a discrepancy in classification based on hits associated with a contig. Defaults to 0.75 (75%). 	|
| --busco_file 	| busco_file 	| Overrides specific organism and taxonomy parameters (next two entries below) in favor of a tab-separated file containing each organism/group of interest and the taxonomic level of the query. 	|
| --organisms 	| organisms 	| A list of organisms/groups to test the BUSCO completeness of matching contigs for. 	|
| --taxonomy_organisms 	| taxonomy_organisms 	| The taxonomic level of the groupings indicated in the list of `--organisms`; also a list. 	|
| --individual_or_summary 	| individual_or_summary 	| Defaults to summary. Whether BUSCO assessment should just be performed for the top organism matches, or whether the list of organisms + their taxonomies or BUSCO file (above parameters) should be used (individual). 	|
| --busco_threshold 	| busco_threshold 	| The threshold for BUSCO completeness for a set of contigs to be considered reasonably BUSCO-complete. 	|
| --create_tax_table 	| create_tax_table 	| If included in a command line argument or set to 1 in a configuration file, the taxonomy table and protein dictionary are produced automatically from the provided reference FASTA file. Defaults to False, then a check is performed. 	|
| --original_tax_table 	| original_tax_table 	| If a new tax table is to be created, the taxonomy table originally included with the database needs to be specified. If the database is downloaded from scratch, this variable is set automatically. Required if `--create-tax-table` is set to True. 	|
| --reformat_tax 	| reformat_tax 	| If included in a command line argument or set to 1 in a configuration file, this means that the original taxonomy table we provided has semicolon-separated taxonomy entries, rather than separate columns for each taxonomic group. Defaults to False, then a check is performed. 	|
| --strain_col_id 	| strain_col_id 	| The column of the provided original taxonomy file which indicates strain name. Defaults to "strain_name". 	|
| --taxonomy_col_id 	| taxonomy_col_id 	| The column of the provided original taxonomy file which indicates strain taxonomy. Defaults to "taxonomy". 	|
| --delimiter 	| delimiter 	| The separator for the taxonomy file entries. Defaults to "/". 	|
| --tax_table 	| tax_table 	| The name of the formatted taxonomy table; defaults to "tax-table.txt.". If this file is not found, it will be generated from the reference FASTA and original taxonomy file. 	|
| --protein_map 	| protein_map 	| The name of the JSON file containing protein correspondences; defaults to "protein-map.json". If this file is not found, it will be generated from the reference FASTA and original taxonomy file. 	|
|  	|  	|  	|

### Example usage

### Community guidelines 

#### How to contribute to `EUKulele`


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
