# Executing EUKulele

SOP for MAG and metatranscriptome analysis with EUKulele. Tutorial details EUKulele usage with protistan metatranscriptome data, explores EUKulele database options, and how to compile output files in R for downstream data analysis.

## (1) Installation & working environment for EUKulele
```
## include pip install or conda install instructions here ##

conda activate EUKulele

EUKulele --help # Pull up help menu
```

## (2) Set up metatranscriptome data & databases

### _(2.1) Data structure_
These data were assembled using megahit, which generates a directory for each assembly with contigs in a file called ```final.contigs.fa```. To run EUKulele and downstream steps, I needed these non-descript file names to have actual sample names. So I moved them to a single assemly directory and renamed them. Additionally, EUKulele requires sample-specific naming for input samples.
_Cite original data from NB_
```
# Location of all assembly files with file suffix ".fa"
assemblies/
```

### _(2.2) Download and view contents of EUKulele databases_

EUKulele defaults to downloading the selected database, unless the user specifies a custom input database. When no database is specified, it will automatically download and use MMETSP. To only download the database, run the following command:
```
# Within the EUKulele conda environment
EUKulele download --reference_dir phylodb
```
Options for ```--referenece_dir``` include: ```mmetsp```, ```phylodb```, ```eukprot```, or ```eukzoo```.

#
# Insert figures of different database compositions.
#

> After downloading the databases, you may want to compare the outputs from all three. And for metatranscriptomics blastx is recommended.


## (3) Running EUKulele
The default usage of EUKulele will look for your assembly files and estimate taxonomy using the MMETSP database.
> Run with defaults, will download MMETSP database in your working directory and run on all files with ```.fa``` ending.
```
EUKulele -m mets -s assemblies/ --n_ext .fa
```

### _(3.1) Additional ways to execute EUKulele_

```
# Run with mmetsp
EUKulele -m mets \
        -s assemblies/ \                                                                                                                   
        --alignment_choice BLAST \
        --reference_dir /vortexfs1/omics/alexander/shu/db/mmetsp/ \
        -o /vortexfs1/scratch/sarahhu/eukulele-output \
        --n_ext .fa

EUKulele -m mets \
        -s assemblies/ \
        --alignment_choice blast \
        --reference_dir /vortexfs1/omics/alexander/shu/db/phylodb/ \
        -o /vortexfs1/scratch/sarahhu/eukulele-output-phylodb \
        --n_ext .fa

EUKulele -m mets \
        -s assemblies/ \
        --alignment_choice blast \
        --reference_dir /vortexfs1/omics/alexander/shu/db/eukprot/ \
        -o /vortexfs1/scratch/sarahhu/eukulele-output-eukprot \
        --n_ext .fa

```

> Same execution, but using a config file
```
EUKulele --config config.yaml
```
The contents of the config file:

```
jobname: run-EUKulele
mets_or_mags: mets
nucleotide_extension: .fa
output: /vortexfs1/scratch/sarahhu/test-eukulele-mmetsp
samples: /vortexfs1/omics/alexander/shu/axial-seamount/tmp

# Database selection
database: mmetsp
reference: /vortexfs1/omics/alexander/shu/db
```

## (4) Explore outputs from EUKulele using R

Import EUKulele taxonomy assignment files.   

### _(4.1) Import one taxonomy assignment file_
**Recommended approach**

* Insert R code here to import from MMETSP, phylodb, eukprot, and eukzoo(?)
* parse and summarize by sample
* explain with and without salmon count files

### _(4.2) With Salmon counts_


### _(4.xx) Compile outputs from multiple outputs_

* Use R to loop and compile all output files into 1 dataframe
* Summarize outputs



* Write shell script to move tiered files to 1 location and change name to sample ID of some kind
* R notebook to provide information on all counts, compile all output data

- Defaults to downloading mmetsp
- option to download first and store, download all and store?
- force overwriting or updating of ref database too, that would be helpful?



## Tutorial

