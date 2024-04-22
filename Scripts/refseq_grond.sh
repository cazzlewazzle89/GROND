#!/bin/sh
VAR_PARALLEL_DL=50 # number of parallel processes to run when downloading assemblies & annotations
VAR_PARALLEL_EDIT=50 # number of parallel processes to run when editing seqid headers
VAR_THREADS_SEQLENGTH=50 # number of threads to use when calculating sequence lengths with seqkit
VAR_THREADS_VSEARCH=50 # number of threads to use then clustering operon sequences with vsearch
VAR_CLUSTERID_VSEARCH=0.999 # vsearch clustering identity
VAR_DB_TAXONDB="/home/cwwalsh/Databases/TaxonKit/" # path to Taxonkit database
VAR_OUTPUT_DIRECTORY=Database_RefSeq_207 # specify output directory name for databases (path must be relative from Current Working Directory)
VAR_TEMP_DIRECTORY=${PWD}/Grond_RefSeq_207/ # specify temp directory to store intermediary files (path must be relative from Current Working Directory)

VAR_SOURCE_DIRECTORY=${PWD}

# make output directory and move shell to run from there
mkdir ${VAR_TEMP_DIRECTORY}
cd ${VAR_TEMP_DIRECTORY}

# make directories to store genomes used for building the database
mkdir -p Genomes_Complete/
mkdir -p Genomes_Incomplete/

# get information for the lastest version of complete prokaryotic genome assemblies from RefSeq
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt \
    -O bacteria_assembly_summary.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt \
    -O archaea_assembly_summary.txt

cat archaea_assembly_summary.txt bacteria_assembly_summary.txt > prokaryote_assembly_summary.txt
rm -f archaea_assembly_summary.txt bacteria_assembly_summary.txt
awk -F '\t' '$11=="latest"' \
    prokaryote_assembly_summary.txt > prokaryote_assembly_summary_latest.txt
awk -F '\t' '$12=="Complete Genome"' \
    prokaryote_assembly_summary_latest.txt > prokaryote_assembly_summary_latestcomplete.txt
awk -F '\t' '$12!="Complete Genome"' \
    prokaryote_assembly_summary_latest.txt > prokaryote_assembly_summary_latestincomplete.txt

# extract ftp directory paths for each assembly
cut -f 20 prokaryote_assembly_summary_latestcomplete.txt > ftpdirpathscomplete
cut -f 20 prokaryote_assembly_summary_latestincomplete.txt > ftpdirpathsincomplete

# compress assembly information
# retaining to extract taxids and keep information about database version
gzip prokaryote_assembly_summary.txt 
gzip prokaryote_assembly_summary_latest.txt 
gzip prokaryote_assembly_summary_latestcomplete.txt 
gzip prokaryote_assembly_summary_latestincomplete.txt

# download nucleotide sequences
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Genomes_Complete/ "ftpdir,file}' \
    ftpdirpathscomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Genomes_Incomplete/ "ftpdir,file}' \
    ftpdirpathsincomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

# simplify assembly filenames
find Genomes_Complete/ -type f -name *.fna.gz | awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' | bash
find Genomes_Incomplete/ -type f -name *.fna.gz | awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' | bash

# edit fasta headers to contain the RefSeq assembly accession
for i in $(find Genomes_Complete/ -type f -name *.fna.gz | sed 's/Genomes_Complete\///' | sed 's/.fna.gz//')
do
    echo mv Genomes_Complete/${i}.fna.gz Genomes_Complete/temp_${i}.fna.gz \; \
        zcat Genomes_Complete/temp_${i}.fna.gz \| \
        sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
        gzip \> Genomes_Complete/${i}.fna.gz \; \
        rm -f Genomes_Complete/temp_${i}.fna.gz
done > do_edit.sh
parallel -a do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f do_edit.sh

for i in $(find Genomes_Incomplete/ -type f -name *.fna.gz | sed 's/Genomes_Incomplete\///' | sed 's/.fna.gz//')
do
    echo mv Genomes_Incomplete/${i}.fna.gz Genomes_Incomplete/temp_${i}.fna.gz \; \
    zcat Genomes_Incomplete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> Genomes_Incomplete/${i}.fna.gz \; \
    rm -f Genomes_Incomplete/temp_${i}.fna.gz
done > do_edit.sh
parallel -a do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f do_edit.sh

# make directories to store annotation information
mkdir Annotation_Complete/
mkdir Annotation_Incomplete/

# download annotation
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Annotation_Complete/ "ftpdir,file}' \
    ftpdirpathscomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Annotation_Incomplete/ "ftpdir,file}' \
    ftpdirpathsincomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

rm -f ftpdirpathscomplete
rm -f ftpdirpathsincomplete

# simplify annotation filenames
find Annotation_Complete/ -type f -name *.gff.gz | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' | bash
find Annotation_Incomplete/ -type f -name *.gff.gz | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' | bash

# extract rRNA genes and edit annotation seqids to contain RefSeq accession and match fasta headers
for i in $(find Annotation_Complete/ -type f -name '*.gff.gz' | sed 's/Annotation_Complete\/// ; s/\.gff\.gz//')
do
    echo zcat Annotation_Complete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> Annotation_Complete/${i}_rrna.gff
done > getRNA.sh
parallel -a getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f getRNA.sh

for i in $(find Annotation_Incomplete/ -type f -name '*.gff.gz' | sed 's/Annotation_Incomplete\/// ; s/\.gff\.gz//')
do
    echo zcat Annotation_Incomplete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> Annotation_Incomplete/${i}_rrna.gff
done > getRNA.sh
parallel -a getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f getRNA.sh

# make directories to store outputs
mkdir Outputs_Complete/
mkdir Outputs_Incomplete/

# make single multifasta files for all assemblies and annotations
find Genomes_Complete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > Outputs_Complete/combined.fna
find Genomes_Incomplete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > Outputs_Incomplete/combined.fna
find Annotation_Complete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > Outputs_Complete/combined_rrna.gff
find Annotation_Incomplete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > Outputs_Incomplete/combined_rrna.gff

# calculate the length of all sequences in mulitfasta
# this info will be used later to discard operons which cross the start/end boundary of the chromosome
seqkit fx2tab Outputs_Complete/combined.fna -n -l -j ${VAR_THREADS_SEQLENGTH} > Outputs_Complete/seq_length.tsv
seqkit fx2tab Outputs_Incomplete/combined.fna -n -l -j ${VAR_THREADS_SEQLENGTH} > Outputs_Incomplete/seq_length.tsv

# run custom R script to filter and reformat annotation
# takes GFF annotation as input
# for every 16S gene, identifies the closest 23S gene with the same seqid and strand values
# identifies unlined rRNA operons and those that cross the chromosome start/end boundary, removes these and writes their details to 'operons_unlinked_or_crossingboundary.tsv'
# writes coordinates of ITS regions to 'combined_ITS.gff'
# writes coordinates of 16S-ITS-23S regions full_operons.gff' for extraction later with bedtools
# writes 'operon_identifiers.txt' with the unique operon identifiers that are substituted into the fasta headers later
Rscript ${VAR_SOURCE_DIRECTORY}/Scripts/refseq_complete_filter_make_gff.R
Rscript ${VAR_SOURCE_DIRECTORY}/Scripts/refseq_incomplete_filter_make_gff.R

# create a multifasta file containing all rRNA operons
bedtools getfasta \
    -fi Outputs_Complete/combined.fna \
    -bed Outputs_Complete/full_operons.gff \
    -fo Outputs_Complete/full_operons.fna \
    -s
rm -f Outputs_Complete/combined.fna.fai

bedtools getfasta \
    -fi Outputs_Incomplete/combined.fna \
    -bed Outputs_Incomplete/full_operons.gff \
    -fo Outputs_Incomplete/full_operons.fna \
    -s
rm -f Outputs_Incomplete/combined.fna.fai

# remove combined assembly and annotation files to save space
# temporarily keeping individual assembly files to do pairwise ANI
rm -f Outputs_Complete/combined.fna Outputs_Incomplete/combined.fna
rm -f Outputs_Complete/combined_rrna.gff Outputs_Combined/combined_rrna.gff

# renames sequences in operon fasta file by giving each operon a unique identifier
cut -f 2 Outputs_Complete/operon_identifiers.txt | \
    sed 's/^/>/' > Outputs_Complete/operon_identifiers_headers.txt
cut -f 2 Outputs_Incomplete/operon_identifiers.txt | \
    sed 's/^/>/' > Outputs_Incomplete/operon_identifiers_headers.txt

# replaces the fasta headers with unique operon identifiers
sed -i 's/>.* />/' Outputs_Complete/operon_identifiers_headers.txt
awk 'NR%2==0' Outputs_Complete/full_operons.fna | \
    paste -d '\n' Outputs_Complete/operon_identifiers_headers.txt - > Outputs_Complete/full_operons_identifiers.fna

sed -i 's/>.* />/' Outputs_Incomplete/operon_identifiers_headers.txt
awk 'NR%2==0' Outputs_Incomplete/full_operons.fna | \
    paste -d '\n' Outputs_Incomplete/operon_identifiers_headers.txt - > Outputs_Incomplete/full_operons_identifiers.fna

# create list of assembly_accessions and taxonIDs from summary file
zcat prokaryote_assembly_summary_latest.txt.gz | \
    cut -f 1,6 > assemblyaccession_taxid.txt

# modify first column to remove the .[0-9] at the end of accession to match prefix added to fasta header and seqid
sed -i 's/\.[0-9]//' assemblyaccession_taxid.txt

# assigning taxonomy lineage to taxids
# output is a tab delim file with one operon per row
taxonkit lineage \
    --data-dir ${VAR_DB_TAXONDB} assemblyaccession_taxid.txt \
    -i 2 | \
    awk '$2!=""' > assemblyaccession_taxid_lineage.txt

# converting operon lineage info to MPA-style
taxonkit reformat \
    -a \
    -i 3 \
    -d ";" \
    -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" \
    -P assemblyaccession_taxid_lineage.txt \
    --data-dir ${VAR_DB_TAXONDB} | \
    cut -f 1,4 > assemblyaccession_taxid_lineage_mpa.txt

# tidying up intermediary files
rm -f assemblyaccession_taxid.txt assemblyaccession_taxid_lineage.txt
rm -f Outputs_Complete/operon_identifiers_headers.txt  Outputs_Complete/full_operons.fna
rm -f Outputs_Incomplete/operon_identifiers_headers.txt Outputs_Incomplete/full_operons.fna

# sorts sequences by length (important for NR clustering, will retain longest sequence in each cluster)
sortbyname.sh \
    in=Outputs_Complete/full_operons_identifiers.fna \
    out=Outputs_Complete/full_operons_sorted.fna \
    length \
    descending

sortbyname.sh \
    in=Outputs_Incomplete/full_operons_identifiers.fna \
    out=Outputs_Incomplete/full_operons_sorted.fna \
    length \
    descending

# combines all operons to make a comprehensive multifasta
mkdir Outputs_Combined/

cat Outputs_Complete/full_operons_identifiers.fna \
    Outputs_Incomplete/full_operons_identifiers.fna > Outputs_Combined/full_operons_identifiers.fna

rm -f Outputs_Complete/full_operons_identifiers.fna
rm -f Outputs_Incomplete/full_operons_identifiers.fna

# removes redundant sequences from Complete Assemblies dataset based on specified identity cutoff
# also extracts relevant cluster information
vsearch \
    --cluster_smallmem Outputs_Complete/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids Outputs_Complete/nrRep.fna \
    --uc Outputs_Complete/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log Outputs_Complete/vsearchlog.txt \
    --consout Outputs_Complete/nrCon.fna

cat Outputs_Complete/clusters.tsv | \
    awk '$1=="S"' > Outputs_Complete/vsearch_centroids.tsv
cat Outputs_Complete/clusters.tsv | \
    awk '$1=="C"' > Outputs_Complete/vsearch_clusters.tsv
cat Outputs_Complete/clusters.tsv | \
    awk '$1=="H"' > Outputs_Complete/vsearch_hits.tsv

# adds operon sequences from Incomplete Assemblies to the NR dataset generated from complete assemblies, reruns clustering, and extracts relevant cluster information
# only adds NR sequences from Incomplete genome assemblies if they are not already represented in the Complete dataset
# also combines the taxonomy tables of these incomplete operons
cat Outputs_Complete/nrRep.fna \
    Outputs_Incomplete/full_operons_sorted.fna > Outputs_Combined/full_operons_sorted_vsearchinput.fna

vsearch \
    --cluster_smallmem Outputs_Combined/full_operons_sorted_vsearchinput.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids Outputs_Combined/nrRep.fna \
    --uc Outputs_Combined/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log Outputs_Combined/vsearchlog.txt \
    --consout Outputs_Combined/nrCon.fna \
    --usersort
cat Outputs_Combined/clusters.tsv | \
    awk '$1=="S"' > Outputs_Combined/vsearch_centroids.tsv
cat Outputs_Combined/clusters.tsv | \
    awk '$1=="C"' > Outputs_Combined/vsearch_clusters.tsv
cat Outputs_Combined/clusters.tsv | \
    awk '$1=="H"' > Outputs_Combined/vsearch_hits.tsv

rm -f Outputs_Combined/full_operons_sorted_vsearchinput.fna

# edit fasta headers of nrCon sequence files
sed -i 's/centroid=// ; s/;seqs=.*//' Outputs_Complete/nrCon.fna 
sed -i 's/centroid=// ; s/;seqs=.*//' Outputs_Combined/nrCon.fna

# move database files to single directory
mkdir -p Database/Full/
mkdir -p Database/NR/

mv Outputs_Combined/full_operons_identifiers.fna Database/Full/full.fna
gzip Database/Full/full.fna

sortbyname.sh \
    in=Outputs_Combined/nrRep.fna \
    out=Database/NR/nrRep.fna.gz

sortbyname.sh \
    in=Outputs_Combined/nrCon.fna \
    out=Database/NR/nrCon.fna.gz

# remove database sequence files from output directories
rm -f Outputs_Combined/nrRep.fna
rm -f Outputs_Combined/nrCon.fna

# generate database taxonomy files and move to database directory
Rscript ${VAR_SOURCE_DIRECTORY}/Scripts/refseq_get_taxonomy.R
mv Outputs_Combined/taxFull.tsv Database/Full/taxFull.tsv
mv Outputs_Combined/taxRep.tsv Database/NR/taxRep.tsv
mv Outputs_Combined/taxLCA.tsv Database/NR/taxLCA.tsv
mv Outputs_Combined/taxMaj.tsv Database/NR/taxMaj.tsv

# calculate mean and median genome size for each taxon in the database (at all 7 taxonomic ranks)
# doing this using either all genomes or complete genomes only
Rscript ${VAR_SOURCE_DIRECTORY}/Scripts/refseq_get_genome_length.R
cp stats_genomelength.tsv Database/

# calculate mean and median rrn operon copy number for each taxon in the database (at all 7 taxonomic ranks)
Rscript ${VAR_SOURCE_DIRECTORY}/Scripts/refseq_get_copy_number.R
cp stats_copynumber.tsv Database/

# move database to current working directory
mv Database/ ${VAR_SOURCE_DIRECTORY}/${VAR_OUTPUT_DIRECTORY}

