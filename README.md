# GROND: Genome-derived Ribosomal OperoN Database
### A quality-checked and publicly available database of full-length 16S-ITS-23S rRNA operon sequences
#### (Formerly FANGORN)
This repository makes available the scripts used to download or build the GROND databases described in the [manuscipt](https://doi.org/10.1099/mgen.0.001255) and available for download [here](https://zenodo.org/records/10889037).  

For a more in-depth comparison of 16S-ITS-23S analysis methods, see our [preprint](https://www.researchsquare.com/article/rs-4006805/v1) (currently under review).  

Please get in touch if you have any comments, issues, or suggestions for improvements.

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
The simplest method is to use the `download` button on the website.  
For downloading with `wget` or `curl` - you can copy the URL of each file, minus the `?download=1` suffix.  
For the sake of convenience, download scripts have been included in this repository - they just contain the `curl` commands ready to go and can be used by simply running `sh download_gtdb_nr.sh`.  

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

`refseq_grond.sh` and `gtdb_grond.sh` are used to build databases using RefSeq and GTDB taxonomy respectively.  
By default they will use:  
- 50 parallel processes to download assemblies and annotations and edit seqid headers  
- 50 threads for VSEARCH clustering  
- 20 threads for identifying rRNA genes in genomes without accompanying annotation using Barrnap  

These settings, along with VSEARCH clustering identity (default: 99.9%) and output directory (default: `./Grond_GTDB/` or `./Grond_RefSeq/`) can be changed by editing the variables at the beginning of the relevant script.  

Both scripts mentioned above will do the majority of the work but will call R scripts to identify _rrn_ operons, perform quality checking (the details of which can be found in the GROND publication), and join _rrn_ operon identifiers to the source genome taxonomy.
R scripts are also used by the GTDB database builder to extract assembly information and identify genomes which are missing annotation.  

Also included in the `Scripts/` folder are scripts to:  
* collate the files necessary to generate the manuscript figures and statistics  
* calculate average nucleotide identity (ANI) between all pairwise combinations of complete genomes used to construct both GTDB and RefSeq databases  
* train a Qiime2 Naive Bayes feature classifier to assign taxonomy to amplicons (generated using all tested primer pairs) using each database - I will make these available shortly  
* calculate intragenomic diversity of complete GTDB genomes  

## Note on Database Construction 

GROND is envisaged as a tool to aid standardisation of 16S-ITS-23S rRNA analysis and allow comparison of results and, as such, building your own version would defeat the purpose.  

That said, if you want to build your own version using the NCBI taxonomy system, make sure you have the most up-to-date version of the taxonomy database. I do this using the commands described in the [TaxonKit manual](https://bioinf.shenwei.me/taxonkit/usage/#before-use).  

This does not download automatically as I often encounter an EOF error while extracting the `taxdump` tarball so it is less hassle to do it manually using the accomanying script `get_taxonkit_DB.sh` before running `refseq_grond.sh`.  
You do not need to do this when using the GTDB taxonomy system as the information is available on [the website](https://gtdb.ecogenomic.org/downloads) and will be automatically downloaded when running `gtdb_grond.sh`.

## Note on Naming

GROND is named after the [wolf-shaped battering ram](https://lotr.fandom.com/wiki/Grond_(battering_ram)) that broke the Great Gate of Minas Tirith during the Battle of the Pelennor Fields, itself named in homage to the [great hammer of Angband](https://lotr.fandom.com/wiki/Grond_(hammer)) wielded by Morgoth, the first Dark Lord.  

It was originally named FANGORN (Full-length Amplicons for the Next Generation Of rRNa analysis) but was renamed to avoid confusion with the much more popular [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html).  

There may still be historic references to FANGORN in the scripts - please reach out if you find these and I will fix them  
