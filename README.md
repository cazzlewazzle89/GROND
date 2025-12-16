# GROND: Genome-derived Ribosomal OperoN Database
### A quality-checked and publicly available database of full-length 16S-ITS-23S rRNA operon sequences
This repository makes available the scripts used to build the GROND databases described in the [manuscipt](https://doi.org/10.1099/mgen.0.001255) and available for download [here](https://zenodo.org/records/17704151).  

For a more in-depth comparison of 16S-ITS-23S analysis methods, see our [publication](https://www.nature.com/articles/s41598-024-83410-7) in Scientific Reports.  

Please get in touch if you have any comments, issues, or suggestions for improvements.

Note: this is an updated version of the database whose construction differs from that described in the original manuscript. 
For release 226 and all future releases - GROND will include only genomes sourced from the Genome Taxonomy Database (GTDB).
For this same set of genomes, two parallel taxonomy annotations are now provided: one based on GTDB lineages and one based on NCBI lineages. This allows users to work within either taxonomic framework without needing to maintain separate genome collections or perform additional remapping steps.

In addition, the dereplication strategy has been expanded to include multiple nucleotide identity thresholds, providing greater flexibility for downstream analyses by allowing users to select genome sets at different levels of redundancy, depending on compute resources and desired phylogenetic resolution.

I plan to update the database in line with each major GTDB release.  
If I am behind the times and don't have a database for the latest GTDB release then let me know by posting an issue.  

## Quick Links
[Database Download](https://github.com/cazzlewazzle89/GROND#database-download)  
[Dependencies](https://github.com/cazzlewazzle89/GROND#dependencies)  
[Usage](https://github.com/cazzlewazzle89/GROND#usage)  
[Note on Database Construction](https://github.com/cazzlewazzle89/GROND#note-on-database-construction)  
[Note on Naming](https://github.com/cazzlewazzle89/GROND#note-on-naming)  

## Database Download

The database is hosted on Zenodo so downloading each file individually is relatively simple.  
The simplest method is download with `wget` or `curl` - you can copy the URL of each file, minus the `?download=1` suffix.  
For example, you would use `wget https://zenodo.org/records/17704151/files/full.fna.gz` to download the full database. 

## Dependencies for Database Construction
Make sure these are in your $PATH

### Command Line Tools
| Software  | Version Tested |
| --- | --- |
| [Barrnap](https://github.com/tseemann/barrnap) | 0.9 |
| [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) | 38.90  |
| [BEDTools](https://github.com/arq5x/bedtools2) | 2.30.0  |
| [csvtk](https://github.com/shenwei356/csvtk) | 0.24.0 |
| [R](https://www.r-project.org/) | 4.1.0  |
| [SeqKit](https://github.com/shenwei356/seqkit) | 2.2.0 |
| [TaxonKit](https://bioinf.shenwei.me/taxonkit/) | 0.8.0  |
| [VSEARCH](https://github.com/torognes/vsearch) | 2.17.1  |

### R Packages

| Package | Version Tested |
| --------|----------------|
| [Tidyverse](https://www.tidyverse.org/) | 1.3.1 |
| [Ape](https://cran.r-project.org/web/packages/ape/index.html) | 5.0 |

The conda environment that I used to build the database can be created using `conda env create -f grond.yml` (the yml is provded in this repository).  
This will create an environment called `grond` which can be loaded with `conda activate grond`.  
This environment contains everything except `csvtk` and `Seqkit`. These are available through bioconda but I used preexisting non-conda installations.

## Usage

`grond.sh` is used to build the database from the latest GTDB release.
By default it will use:  
- 50 parallel processes to download assemblies and annotations and edit seqid headers  
- 50 threads for VSEARCH clustering  
- 20 threads for identifying rRNA genes in genomes using Barrnap  

These settings, along with output directory can be changed by editing the variables at the beginning of the relevant script.  

Update from the original version: all dereplicated databases will be constructed by running `grond.sh` - I removed the option to specify the dereplication threshold. 

This script will do the majority of the work but will call R scripts to identify _rrn_ operons, perform quality checking (the details of which can be found in the GROND publication), and join _rrn_ operon identifiers to the source genome taxonomy.
R scripts are also used by the GTDB database builder to extract assembly information.

## Note on Database Construction 

GROND is envisaged as a tool to aid standardisation of 16S-ITS-23S rRNA analysis and allow comparison of results and, as such, building your own version would defeat the purpose.  

## Note on Naming

GROND is named after the [wolf-shaped battering ram](https://lotr.fandom.com/wiki/Grond_(battering_ram)) that broke the Great Gate of Minas Tirith during the Battle of the Pelennor Fields, itself named in homage to the [great hammer of Angband](https://lotr.fandom.com/wiki/Grond_(hammer)) wielded by Morgoth, the first Dark Lord.  

It was originally named FANGORN (Full-length Amplicons for the Next Generation Of rRNa analysis) but was renamed to avoid confusion with the much more popular [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html).  

There may still be historic references to FANGORN in the scripts - please reach out if you find these and I will fix them  
