# GROND: Genome-derived Ribosomal OperoN Database
![Zenodo Downloads](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fzenodo.org%2Fapi%2Frecords%2F17704151&label=downloads&query=%24.stats.downloads&color=blue) 
### A quality-checked and publicly available database of full-length 16S-ITS-23S rRNA operon sequences
This repository makes available the scripts used to build the GROND databases described in the [manuscript](https://doi.org/10.1099/mgen.0.001255) and available for download [here](https://zenodo.org/records/21278902).

For a more in-depth comparison of 16S-ITS-23S analysis methods, see our [publication](https://www.nature.com/articles/s41598-024-83410-7) in Scientific Reports.

Please get in touch if you have any comments, issues, or suggestions for improvements.

> **Note:** This is an updated version of the pipeline whose construction differs from that described in the original manuscript.
> The pipeline now accepts a user-supplied manifest file rather than downloading genomes directly from GTDB/NCBI. This makes it easier to build custom databases from any set of genome assemblies with pre-computed barrnap annotations.
>
> In addition to the 16S-ITS-23S rRNA operon databases, companion **16S rRNA** and **23S rRNA** sequence databases are now generated at independent dereplication thresholds. All databases are dereplicated — there is no longer an unduplicated "full" database; instead, the least-restrictive threshold (100% identity) serves as the base database.
>
> Two parallel taxonomy annotations continue to be provided: one based on GTDB lineages and one based on NCBI lineages.

I plan to update the database in line with each major GTDB release.
If I am behind the times and don't have a database for the latest GTDB release then let me know by posting an issue.

## Quick Links
- [Database Download](#database-download)  
- [Database Contents](#database-contents)  
- [Dependencies](#dependencies-for-database-construction)  
- [Preparing the Manifest](#preparing-the-manifest)  
- [Usage](#usage)  
- [Note on Database Construction](#note-on-database-construction)  
- [Note on Naming](#note-on-naming)  

## Database Download

The database is hosted on Zenodo so downloading each file individually is relatively simple.
The simplest method is to use `wget` or `curl` to target the URL of each file, minus the `?download=1` suffix.
For example, you would use `wget https://zenodo.org/records/17704151/files/operon_nrRep_1.0.fna.gz` to download the 100% identity dereplicated operon database.

## Database Contents

All databases are named by sequence type and percent-identity (PID) threshold used for dereplication.

### rRNA Operon (16S-ITS-23S) Databases

| File | PID Threshold | Description | Sequence Count |
| --- | --- | --- | --- |
| `operon_nrRep_1.0.fna.gz` | 100% | Representative sequences - all unique operons | 361,183 |
| `operon_nrRep_0.999.fna.gz` | 99.9% | Representative sequences | 240,583 |
| `operon_nrRep_0.99.fna.gz` | 99% | Representative sequences | 90,986 |
| `operon_nrRep_0.97.fna.gz` | 97% | Representative sequences | 53,076 |
| `operon_nrRep_0.95.fna.gz` | 95% | Representative sequences | 36,688 |
| `operon_nrRep_0.90.fna.gz` | 90% | Representative sequences | 17,547 |
| `operon_nrCon_*.fna.gz` | (same thresholds) | Consensus sequences for each cluster | * |

### 16S rRNA Databases

| File | PID Threshold | Description | Sequence Count |
| --- | --- | --- | --- |
| `16S_nrRep_1.0.fna.gz` | 100% | Representative sequences - all unique genes | 174,303 |
| `16S_nrRep_0.99.fna.gz` | 99% | Representative sequences | 46,449 |
| `16S_nrRep_0.98.fna.gz` | 98% | Representative sequences | 32,035 |
| `16S_nrRep_0.97.fna.gz` | 97% | Representative sequences | 24,254 |
| `16S_nrCon_*.fna.gz` | (same thresholds) | Consensus sequences for each cluster | * |

### 23S rRNA Databases

| File | PID Threshold | Description | Sequence Count |
| --- | --- | --- | --- |
| `23S_nrRep_1.0.fna.gz` | 100% | Representative sequences - all unique genes | 253,485 |
| `23S_nrRep_0.99.fna.gz` | 99% | Representative sequences | 56,680 |
| `23S_nrRep_0.98.fna.gz` | 98% | Representative sequences | 41,016 |
| `23S_nrRep_0.97.fna.gz` | 97% | Representative sequences | 32,347 |
| `23S_nrCon_*.fna.gz` | (same thresholds) | Consensus sequences for each cluster | * |

### Taxonomy Files

For each dereplication threshold and each sequence type, three taxonomy assignment schemes are provided for both GTDB and NCBI frameworks:

| Scheme | Description |
| --- | --- |
| `taxRep` | Taxonomy of the cluster representative sequence |
| `taxLCA` | Lowest Common Ancestor of all sequences in the cluster |
| `taxMaj` | Lowest rank with a single majority agreement across all sequences in the cluster |

A `taxFull` file is also provided mapping every sequence in the source dataset to its source genome taxonomy.

### Cluster Membership Files

| File | Description |
| --- | --- |
| `operon_clusterMap.tsv`, `16S_clusterMap.tsv`, `23S_clusterMap.tsv` | Every sequence in the full (pre-dereplication) dataset, whether it came from a complete or incomplete genome (`Source`), and the OperonID of its representative sequence in the `nrRep` database at each threshold (`Rep_<pid>` columns) |

### Supplementary Statistics

| File | Description |
| --- | --- |
| `master_rrna.tsv` | Full coordinate table for every identified rRNA operon — includes seqid, OperonID, strand, and positions of the 16S gene, ITS region, and 23S gene |
| `stats_genomelength_gtdb.tsv` | Mean/median genome size per taxon (GTDB) |
| `stats_genomelength_ncbi.tsv` | Mean/median genome size per taxon (NCBI) |
| `stats_copynumber_gtdb.tsv` | Mean/median rrn copy number per taxon (GTDB) |
| `stats_copynumber_ncbi.tsv` | Mean/median rrn copy number per taxon (NCBI) |

## Dependencies for Database Construction

The software below must be in your `$PATH`. The simplest way to install everything is with the supplied conda environment.

### Installation

Clone the repository and create the `grond` environment from `grond.yml` (`mamba` can be substituted for `conda` for a faster install):

```bash
git clone https://github.com/cazzlewazzle89/GROND.git
cd GROND
conda env create -f grond.yml
conda activate grond
```

Alternatively, the environment can be created without the YAML file:

```bash
conda create -n grond -c conda-forge -c bioconda \
    "python>=3.9" "pandas>=1.3" "numpy>=1.21" \
    "barrnap>=0.9" "bedtools>=2.30" "bbmap>=38.90" "seqkit>=2.2" "vsearch>=2.17"
```

Conda also installs HMMER (required by Barrnap) and Java (required by BBTools) as dependencies. The shell scripts additionally use standard Unix utilities (`bash`, `awk`, `sed`, `gunzip`, `xargs`), which are available on any Linux system. Linux is recommended; some Bioconda packages may not be available for Apple Silicon Macs.

### Command Line Tools

| Software | Version Tested | Notes |
| --- | --- | --- |
| [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) (conda: `bbmap`) | 38.90 | Used for `sortbyname.sh` |
| [BEDTools](https://github.com/arq5x/bedtools2) | 2.30.0 | Used for sequence extraction |
| [Barrnap](https://github.com/tseemann/barrnap) | 0.9 | Run externally; output GFF paths supplied via manifest |
| [SeqKit](https://github.com/shenwei356/seqkit) | 2.2.0 | Used for sequence length calculation |
| [VSEARCH](https://github.com/torognes/vsearch) | 2.17.1 | Used for dereplication clustering |

### Python Requirements

Python 3.9+ with the following packages:

| Package | Notes |
| --- | --- |
| [pandas](https://pandas.pydata.org/) | Core data manipulation |
| [numpy](https://numpy.org/) | Numerical operations |

> **Note:** R and all R package dependencies have been removed. The pipeline now uses Python exclusively for all scripting tasks.

## Preparing the Manifest

The pipeline requires a tab-separated manifest file (`manifest.tsv`) with the following columns:

| Column | Description |
| --- | --- |
| `genomeID` | A unique identifier for the genome (used as a prefix in all output sequence IDs) |
| `path_to_fna` | Path to the genome assembly FASTA file (`.fna` or `.fna.gz`) |
| `path_to_gff` | Path to the Barrnap output GFF file (`.gff` or `.gff.gz`) |
| `gtdb_taxonomy` | GTDB lineage string (semicolon-separated, e.g. `d__Bacteria;p__Bacillota;...`) |
| `ncbi_taxonomy` | NCBI lineage string (semicolon-separated) |
| `genome_completeness` | Either `complete` or `incomplete` |

Example:
```tsv
genomeID	path_to_fna	path_to_gff	gtdb_taxonomy	ncbi_taxonomy	genome_completeness
GCA_000001405	/data/genomes/GCA_000001405.fna.gz	/data/gff/GCA_000001405_rrna.gff	d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus pneumoniae	d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus pneumoniae	complete
```

Both absolute and relative paths are supported for `path_to_fna` and `path_to_gff`. Relative paths are resolved from the directory where `grond.sh` is invoked.

### Filling in Missing GFF Files

If some manifest genomes have a FASTA assembly but no Barrnap GFF, `fill_missing_gff.sh` runs Barrnap on them (kingdom `arc` for GTDB `d__Archaea` genomes and `bac` otherwise), writes the GFFs to a local directory, and writes an updated manifest pointing at them. The input manifest is not modified, rows with a missing FASTA are left unchanged, and existing GFFs in the output directory are reused so interrupted runs can be resumed.

```bash
bash SCRIPTS/fill_missing_gff.sh -i manifest.tsv -o manifest_updated.tsv -d barrnap_gff -t 4 -j 6
bash SCRIPTS/grond.sh manifest_updated.tsv
```

`-t` sets threads per Barrnap job and `-j` the number of parallel jobs; `-T` sets the temporary directory used for decompressing `.gz` assemblies (default `$TMPDIR` or `/tmp`). Genomes where Barrnap fails are listed in `<gff_dir>/failed_genomes.tsv`.

### Sequence ID Format

All output sequences are assigned a unique identifier in the format:

```
{genomeID}__{original_seqid}__{operon_number}
```

Where `operon_number` is a 1-based index incremented per contig, ensuring globally unique IDs across all genomes in the database.

## Usage

```bash
bash SCRIPTS/grond.sh manifest.tsv
```

The script will:
1. Parse the manifest and construct a combined multi-FASTA and combined GFF file with genome-prefixed sequence IDs.
2. Calculate sequence lengths with SeqKit.
3. Run the GFF filtering Python script to identify valid 16S-ITS-23S operons, filter out partial/boundary-crossing/unlinked operons, and extract GFF coordinates for operons, 16S genes, and 23S genes.
4. Extract sequences using BEDTools and rename headers with unique IDs.
5. Perform prioritised two-pass dereplication using VSEARCH:
   - Complete genomes are clustered first to ensure they are preferentially selected as cluster representatives.
   - Incomplete genome sequences are then added and a second clustering pass is run using `--usersort` to maintain priority.
6. Assign taxonomy (taxRep, taxLCA, taxMaj) for each database using both GTDB and NCBI schemes.
7. Calculate per-taxon genome size and rRNA copy number statistics.
8. Collect all output files into the output directory (configurable via `VAR_OUTPUT_DIRECTORY`).

### Configurable Variables

Edit the variables at the top of `grond.sh` to adjust:

| Variable | Default | Description |
| --- | --- | --- |
| `VAR_THREADS_SEQLENGTH` | 24 | Threads for SeqKit sequence length calculation |
| `VAR_THREADS_VSEARCH` | 24 | Threads for VSEARCH clustering |
| `VAR_OUTPUT_DIRECTORY` | `R232` | Output directory name |
| `VAR_TEMP_DIRECTORY` | `${PWD}/TEMP` | Temporary working directory |
| `VAR_ALLOW_MISSING` | 0 | If 0, abort when any manifest genome has a missing FASTA/GFF file (listed in `skipped_genomes.tsv`); set to 1 to continue without them |

## Scripts

| Script | Language | Description |
| --- | --- | --- |
| `grond.sh` | Bash | Main pipeline orchestration script |
| `filter_make_gff.py` | Python | GFF filtering, operon identification, and coordinate extraction |
| `get_taxonomy.py` | Python | Taxonomy assignment (taxRep, taxLCA, taxMaj) for all sequence types and thresholds |
| `get_genome_stats.py` | Python | Per-taxon genome size and rRNA copy number statistics |
| `replace_headers.py` | Python | FASTA header replacement with unique operon/gene IDs |
| `fill_missing_gff.sh` | Bash | Runs Barrnap for manifest genomes with a FASTA but no GFF and writes an updated manifest |
| `tidy_largefiles.sh` | Bash | Removes the TEMP working directory left behind by an aborted or interrupted run |

## Note on Database Construction

GROND is envisaged as a tool to aid standardisation of 16S-ITS-23S rRNA analysis and allow comparison of results and, as such, building your own version would defeat the purpose.

However, the pipeline scripts are provided here to allow inspection of the construction methodology and to enable building of custom databases from user-supplied genome sets where that is scientifically appropriate (e.g. for organism-specific or environment-specific studies).

## Note on Naming

GROND is named after the [wolf-shaped battering ram](https://lotr.fandom.com/wiki/Grond_(battering_ram)) that broke the Great Gate of Minas Tirith during the Battle of the Pelennor Fields, itself named in homage to the [great hammer of Angband](https://lotr.fandom.com/wiki/Grond_(hammer)) wielded by Morgoth, the first Dark Lord.

It was originally named FANGORN (Full-length Amplicons for the Next Generation Of rRNa analysis) but was renamed to avoid confusion with the much more popular [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html).

There may still be historic references to FANGORN in the scripts - please reach out if you find these and I will fix them.
