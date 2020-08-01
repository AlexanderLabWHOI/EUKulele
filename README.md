# EUKulele

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/AlexanderLabWHOI/EUKulele.svg?branch=master)](https://travis-ci.com/AlexanderLabWHOI/EUKulele)

![EUKulele-logo](/.infrastructure/eukulele-logo.png)

## Formalizing environmental eukaryotic taxonomic assignment

### About EUKulele
`EUKulele` is a Python program for taxonomic annotation of microbes in metatranscriptomic and metagenomic samples, with special emphasis on eukaryote discovery. `EUKulele` can be downloaded from [PyPI](https://pypi.org/), or it may be downloaded via `conda` and used as a command-line program. The software includes three major features:
- Database setup and formatting
- Database creation, alignment, and taxonomic estimation
- Assessment of the BUSCO completeness of subsets of contigs at each taxonomic level

### Prerequisites for running EUKulele
In principle, there are two prerequisites for running the software:
1. Metagenomic or metatranscriptomic sample files (unless using the provided sample data)
2. A database to align the contigs from the metagenome/metatranscriptome to

Three databases are supported by default from within `EUKulele`, and may be downloaded and formatted automatically if the user chooses (or if another reference directory is not specified/does not exist):
- [PhyloDB](https://drive.google.com/drive/u/0/folders/0B-BsLZUMHrDQfldGeDRIUHNZMEREY0g3ekpEZFhrTDlQSjQtbm5heC1QX2V6TUxBeFlOejQ)
- [EukProt](https://figshare.com/articles/EukProt_a_database_of_genome-scale_predicted_proteins_across_the_diversity_of_eukaryotic_life/12417881/2)
- [MMETSP](https://zenodo.org/record/1212585#.Xw3PoJNKhTZ)

### Downloading and configuring EUKulele

When used as a Python package, `EUKulele` will independently attempt to install dependent software for you. However, this is prone to error based on the idiosyncracies of your system. Hence, it is recommended that you create and activate a `conda` environment for `EUKulele`, by running the command 

```
conda env create -f src/EUKulele/EUKulele-env.yaml
```

from the base directory, and then typing `conda activate EUKulele`. When `EUKulele` is used as a `conda` package, all dependent software will be checked for automatically.

### Basic usage

#### Python module

Inside of a Python script, include `import EUKulele` in the header.

Then, you may execute `EUKulele` as a Python package within `Python` using: 

```
EUKulele.eukulele(config=config_file)
```

if you have a configuration file to specify (replace the `config_file` variable with this path), or with:

```
EUKulele.eukulele(string_arguments=string_of_arguments)
```

where `string_of_arguments` is a string containing the non-default `EUKulele` options you wish to specify.

### Full list of EUKulele parameters

Parameters may either be passed as arguments directly to EUKulele, or may be specified in a configuration `YAML` file. If running `EUKulele` as a Python package, you can call the function `eukulele()` with two potential arguments, `config` or `string_arguments`, one of which is required, and then a third argument `command_line`, which should never be set manually and is only relevant if the script has been run from the command line. 

```
eukulele(config="", string_arguments="", )
```

where the `config` argument allows you to list the path of the configuration `YAML` file that you wish to use, while the `string_arguments` option allows you to create a string containing all of the arguments that you would use normally using the command-line options specified below. 

`EUKulele` can be run from the command line using: `EUKulele [arguments]`, once installed. The following arguments are supported, either as command-line arguments, as a string of arguments passed to the `eukulele()` function within the `EUKulele` Python module, or as part of a confirmation `YAML` file, with names specified in the table. A configuration file may also be provided when `EUKulele` is run from the command line, when the `--config` flag is used with a configuration file path.

| Flag 	| Configuration File Entry 	| Meaning 	|
|-	|-	|-	|
| --config 	| N/A 	| The path to the configuration file which should be used to retrieve the equivalent of command-line arguments. 	|
| -m/--mets_or_mags 	| mets_or_mags 	| A required flag to indicate whether metatranscriptomic ("mets") or metagenomic ("mags") samples are being used as input. 	|
| -s/--sample_dir 	| samples 	| A required flag to indicate where the samples (metagenomic or metatranscriptomic, depending on "mets_or_mags" flag) are located. 	|
| -o/--out_dir 	| output 	| The path to the directory where output will be stored. Defaults to a folder called `output` in the present working directory. 	|
| --reference_dir 	| reference 	| A flag to indicate where the reference FASTA is stored, or a keyword argument for the dataset to be downloaded and used. Only used if not downloading automatically. 	|
| --ref_fasta 	| ref_fasta 	| The name of the reference FASTA file in `reference_dir`; defaults to reference.pep.fa if not specified, or is set according to the downloaded file if using a keyword argument. 	|
| --database 	| database 	| An optional additional argument for specifying the database name. 	|
| --nucleotide_extension/--n_ext 	| nucleotide_extension 	| The file extension for samples in nucleotide format (metatranscriptomes). Defaults to .fasta. 	|
| --protein_extension/--p_ext 	| protein_extension 	| The file extension for samples in protein format (metatranscriptomes). Defaults to .faa. 	|
| --ref_fasta_ext 	| ref_fasta_ext 	| The file extension for `ref_fasta`, if applicable. Defaults to .fasta. 	|
| -f/--force_rerun 	| force_rerun 	| If included in a command line argument or set to 1 in a configuration file, this argument forces all steps to be re-run, regardless of whether output is already present. 	|
| --use_salmon_counts 	| use_salmon_counts 	| If included in a command line argument or set to 1 in a configuration file, this argument causes classifications to be made based both on number of classified transcripts and by counts. 	|
| --salmon_dir 	| salmon_dir 	| If `--use_salmon_counts` is true, this must be specified, which is the directory location of the `salmon` output/quantification files. 	|
| --names_to_reads 	| names_to_reads 	| A file that creates a correspondence between each transcript name and the number of `salmon`-quantified reads. Can be generated manually via the `names_to_reads.py` script, or will be generated automatically if it does not exist. 	|
| --transdecoder_orfsize 	| transdecoder_orfsize 	| The minimum cutoff size for an open reading frame (ORF) detected by `TransDecoder`. 	||
| --alignment_choice 	| alignment_choice 	| A choice of aligner to use, currently `BLAST` or `DIAMOND`. 	|
| --cutoff_file 	| cutoff_file 	| A `YAML` file, provided in `src/EUKulele/static/`, that contains the percent identity cutoffs for various taxonomic classifications. Any path may be provided here to a user-specified file. 	|
| --filter_metric 	| filter_metric 	| Either evalue, pid, or bitscore (default evalue) - the metric to be used to filter hits based on their quality prior to taxonomic estimation. 	|
| --consensus_cutoff 	| consensus_cutoff 	| The value to be used to decide whether enough of the taxonomic matches are identical to overlook a discrepancy in classification based on hits associated with a contig. Defaults to 0.75 (75%). 	|
| --busco_file 	| busco_file 	| Overrides specific organism and taxonomy parameters (next two entries below) in favor of a tab-separated file containing each organism/group of interest and the taxonomic level of the query. 	|
| --organisms 	| organisms 	| A list of organisms/groups to test the BUSCO completeness of matching contigs for. 	|
| --taxonomy_organisms 	| taxonomy_organisms 	| The taxonomic level of the groupings indicated in the list of `--organisms`; also a list. 	|
| --individual_or_summary 	| individual_or_summary 	| Defaults to summary. Whether BUSCO assessment should just be performed for the top organism matches, or whether the list of organisms + their taxonomies or BUSCO file (above parameters) should be used (individual). 	|
| --busco_threshold 	| busco_threshold 	| The threshold for BUSCO completeness for a set of contigs to be considered reasonably BUSCO-complete. 	|
| --tax_table 	| tax_table 	| The name of the formatted taxonomy table; defaults to "tax-table.txt.". If this file is not found, it will be generated from the reference FASTA and original taxonomy file. 	|
| --protein_map 	| protein_map 	| The name of the JSON file containing protein correspondences; defaults to "protein-map.json". If this file is not found, it will be generated from the reference FASTA and original taxonomy file. 	|
|  	|  	|  	|

### Output organization 

- here, put an example of what your output directory will look like if you run EUKulele using (a) mags or (b) mets

### Running a sample

To test drive EUKulele, you can run the following:

```
wget https://www.dropbox.com/s/l4kvbpqftdad5ib/sample_eukulele.tar.gz?dl=1
tar xzf sample_eukulele.tar.gz?dl=1
```

Then `cd` into `sample_EUKulele`. A problem with tarring across systems and through Dropbox requires you to run: 

```
rm samples_MAGs/._sample_*
```

Create a conda environment and activate it, using:

```
conda env create -f EUKulele-env.yaml
conda activate EUKulele
```

Then download `EUKulele` via `pip` using:

```
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps EUKulele
```

If any dependency is not satisfied, you can install it manually using `pip install <requirement> --user`. 

Once everything is setup and all dependencies satisfied, from within `sample_EUKulele`, run:

```
EUKulele --config curr_config.yaml
```

Then check the folder `test_out_23July` in the current directory.

### Community guidelines 

#### How to contribute to `EUKulele`
If you are interested in modifying `EUKulele`, you may fork the project for your own use, as detailed in the [MIT License](https://github.com/AlexanderLabWHOI/EUKulele/blob/master/LICENSE) we have adopted for the project. In order to contribute, please contact the developers via Arianna Krinos (akrinos (at) mit (dot) edu) after making the desired changes, after which a pull request may be submitted. 

#### Submitting an issue
If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please submit an issue on the [Issues tab](https://github.com/AlexanderLabWHOI/EUKulele/issues) of the project `GitHub` repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

#### Asking for help
If you have questions about how to use `EUKulele`, or would like to seek out collaborations related to this project, you may contact Arianna Krinos at akrinos (at) mit (dot) edu. 

### Acknowledgments

Authors: Arianna Krinos, Sarah Hu, and Harriet Alexander.

We gratefully acknowledge Natalie Cohen, who contributed significantly to project development and testing. 

## TO DO

- Create a link for downloading the MMETSP formatted database from Zenodo etc.
- "We might consider looking into BIOM files but this might be overkill at this point. Helps integrate with an ecosystem of other tools though (e.g. qiime2)"
- "Amalgamate all the hits for the proteins and assess the LCA based on all hits. % cutoffs? How to do this?"
- Output all proteins with top hits and proteins that do not have good matches. (check in on what this meant)
- Stats on next best hit?