# EUKulele

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/AlexanderLabWHOI/EUKulele.svg?branch=master)](https://travis-ci.com/AlexanderLabWHOI/EUKulele)
[![Coverage Status](https://coveralls.io/repos/github/AlexanderLabWHOI/EUKulele/badge.svg?branch=master)](https://coveralls.io/github/AlexanderLabWHOI/EUKulele?branch=master)
[![Documentation Status](https://readthedocs.org/projects/eukulele/badge/?version=latest)](https://eukulele.readthedocs.io/en/latest/?badge=latest)
[![Read the Docs](https://img.shields.io/badge/read-the%20docs-green)](https://eukulele.readthedocs.io/en/latest/)

<p align="center">
  <img width = "500" src="/.infrastructure/eukulele-logo.png"/>
</p>

## Formalizing environmental eukaryotic taxonomic assignment

### About EUKulele

`EUKulele` is a Python program for taxonomic annotation of microbes in metatranscriptomic and metagenomic samples, with special emphasis on eukaryote discovery. `EUKulele` can be downloaded from [PyPI](https://pypi.org/), or it may be downloaded via `conda` and used as a command-line program. The software includes four major features:
- Database setup and formatting
- Database creation, alignment, and taxonomic estimation
- Assessment of the BUSCO completeness of subsets of contigs at each taxonomic level
- Assessment of taxonomic classification using only BUSCO-identified core eukaryotic genes

### Prerequisites for running EUKulele

In principle, there are two prerequisites for running the software:
1. Metagenomic or metatranscriptomic sample files (unless using the provided sample data)
2. A database to align the contigs from the metagenome/metatranscriptome to

Three databases are supported by default from within `EUKulele`, and may be downloaded and formatted automatically if the user chooses (or if another reference directory is not specified/does not exist):
- [PhyloDB](https://drive.google.com/drive/u/0/folders/0B-BsLZUMHrDQfldGeDRIUHNZMEREY0g3ekpEZFhrTDlQSjQtbm5heC1QX2V6TUxBeFlOejQ)
- [EukProt](https://figshare.com/articles/EukProt_a_database_of_genome-scale_predicted_proteins_across_the_diversity_of_eukaryotic_life/12417881/2)
- [MMETSP](https://zenodo.org/record/1212585#.Xw3PoJNKhTZ)

### Basic usage

If installed either with pip or `conda`, `EUKulele` can be invoked via::

```
EUKulele <arguments>
```
    
Where the minimal command would be

```
EUKulele --mets_or_mags <choice of data type> --sample_dir <where samples are located>
```

See the [documentation](https://eukulele.readthedocs.io/en/latest/) for further details.

### Community guidelines 

#### How to contribute to `EUKulele`
If you are interested in modifying `EUKulele`, you may fork the project for your own use, as detailed in the [MIT License](https://github.com/AlexanderLabWHOI/EUKulele/blob/master/LICENSE) we have adopted for the project. In order to contribute, please contact the developers via Arianna Krinos (akrinos (at) mit (dot) edu) after making the desired changes, after which a pull request may be submitted. 

#### Submitting an issue
If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please submit an issue on the [Issues tab](https://github.com/AlexanderLabWHOI/EUKulele/issues) of the project `GitHub` repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

#### Asking for help
If you have questions about how to use `EUKulele`, or would like to seek out collaborations related to this project, you may contact Arianna Krinos at akrinos (at) mit (dot) edu. 

### Acknowledgments

Authors: Arianna Krinos, Sarah Hu, Natalie Cohen, and Harriet Alexander.
