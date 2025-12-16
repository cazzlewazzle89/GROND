#!/bin/sh

#####
# updated version of GROND
# uses only GTDB genomes (for simplicity)
# identifies all rRNA genes with barrnap (for consistency)
# also uses a few different dereplication identities to create multiple NR versions
# taxonomy comes from GTDB-provided taxonomy (both GTDB and NCBI lineages)
#####

VAR_PARALLEL_DL=50 # number of parallel processes to run when downloading assemblies & annotations
VAR_PARALLEL_EDIT=50 # number of parallel processes to run when editing seqid headers
VAR_THREADS_SEQLENGTH=50 # number of threads to use when calculating sequence lengths with seqkit
VAR_THREADS_VSEARCH=50 # number of threads to use then clustering operon sequences with vsearch
VAR_BARRNAP_THREADS=20
VAR_OUTPUT_DIRECTORY=R226 # specify output directory name for databases (path must be relative from Current Working Directory)
VAR_TEMP_DIRECTORY=${PWD}/TEMP # specify temp directory to store intermediary files (path must be relative from Current Working Directory)

VAR_SOURCE_DIRECTORY=${PWD}

# make output directory and move shell to run from there
mkdir ${VAR_TEMP_DIRECTORY}
cd ${VAR_TEMP_DIRECTORY}

# download archaea and bacteria metadata and combine into a single file
wget -c https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz \
    -O ar53_metadata.tsv.gz \
    --no-check-certificate
gunzip ar53_metadata.tsv.gz

wget -c https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz \
    -O bac120_metadata.tsv.gz \
    --no-check-certificate
gunzip bac120_metadata.tsv.gz

mv ar53_metadata.tsv metadata.tsv
sed -e '1d' bac120_metadata.tsv >> metadata.tsv
rm -f bac120_metadata.tsv

csvtk cut -t -f accession,gtdb_taxonomy,ncbi_taxonomy metadata.tsv > taxonomy.tsv

# download genbank and refseq assembly information
wget -c ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt \
    -O assembly_summary_genbank.txt
wget -c ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt \
    -O assembly_summary_refseq.txt

# get assembly info for each assembly in the GTDB database
grep '^GB_GCA' taxonomy.tsv | \
    cut -f 1 | \
        sed 's/GB_//' > genbank_gtdb_accessions.txt
grep '^RS_GCF' taxonomy.tsv \
    | cut -f 1 | \
        sed 's/RS_//' > refseq_gtdb_accessions.txt
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_assembly_info.R

# extract ftp directory paths for each assembly
cut -f 20 assembly_summary_genbank_gtdb_complete.txt > ftpdirpaths_genbank_complete
cut -f 20 assembly_summary_genbank_gtdb_incomplete.txt > ftpdirpaths_genbank_incomplete
cut -f 20 assembly_summary_refseq_gtdb_complete.txt > ftpdirpaths_refseq_complete
cut -f 20 assembly_summary_refseq_gtdb_incomplete.txt > ftpdirpaths_refseq_incomplete

# compress assembly information
# retaining to keep information about database version
gzip assembly_summary_genbank.txt 
gzip assembly_summary_refseq.txt 
gzip assembly_summary_genbank_gtdb_complete.txt 
gzip assembly_summary_genbank_gtdb_incomplete.txt 
gzip assembly_summary_refseq_gtdb_complete.txt
gzip assembly_summary_refseq_gtdb_incomplete.txt

# remove files listing assession numbers
rm -f genbank_gtdb_accessions.txt refseq_gtdb_accessions.txt

# make directories to store genomes used for building the database
mkdir -p Genomes_Complete/
mkdir -p Genomes_Incomplete/

# download nucleotide sequences
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -c -P Genomes_Complete/ "ftpdir,file}' \
    ftpdirpaths_genbank_complete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -c -P Genomes_Complete/ "ftpdir,file}' \
    ftpdirpaths_refseq_complete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -c -P Genomes_Incomplete/ "ftpdir,file}' \
    ftpdirpaths_genbank_incomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -c -P Genomes_Incomplete/ "ftpdir,file}' \
    ftpdirpaths_refseq_incomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

# simplify assembly filenames
find Genomes_Complete/ -type f -name '*.fna.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' > rename.sh
sh rename.sh

find Genomes_Incomplete/ -type f -name '*.fna.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' > rename.sh
sh rename.sh
rm -f rename.sh

# edit fasta headers to contain the RefSeq or GenBank assembly accession
find Genomes_Complete/ -type f -name '*.fna.gz' | sed 's/Genomes_Complete\///' | sed 's/.fna.gz//' > genomes_complete.list
for i in $(cat genomes_complete.list)
do

    echo mv Genomes_Complete/${i}.fna.gz Genomes_Complete/temp_${i}.fna.gz \; \
        zcat Genomes_Complete/temp_${i}.fna.gz \| \
            sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
                gzip \> Genomes_Complete/${i}.fna.gz \; \
                    rm -f Genomes_Complete/temp_${i}.fna.gz

done > do_edit.sh
parallel -a do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f do_edit.sh 

find Genomes_Incomplete/ -type f -name '*.fna.gz' | sed 's/Genomes_Incomplete\///' | sed 's/.fna.gz//' > genomes_incomplete.list
for i in $(cat genomes_incomplete.list)
do

    echo mv Genomes_Incomplete/${i}.fna.gz Genomes_Incomplete/temp_${i}.fna.gz \; \
        zcat Genomes_Incomplete/temp_${i}.fna.gz \| \
            sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
                gzip \> Genomes_Incomplete/${i}.fna.gz \; \
                    rm -f Genomes_Incomplete/temp_${i}.fna.gz

done > do_edit.sh
parallel -a do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f do_edit.sh 

# identify rRNA genes in genomes
mkdir Annotation_Complete/
mkdir Annotation_Incomplete/

for i in $(cat genomes_complete.list)
do

    zcat Genomes_Complete/${i}.fna.gz | \
        barrnap --threads ${VAR_BARRNAP_THREADS} \
            --reject 0.8 > Annotation_Complete/${i}_rrna.gff

done

for i in $(cat genomes_incomplete.list)
do
    zcat Genomes_Incomplete/${i}.fna.gz | \
        barrnap --threads ${VAR_BARRNAP_THREADS} \
            --reject 0.8 > Annotation_Incomplete/${i}_rrna.gff
done

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

# run R script to filter and reformat annotation
# takes GFF annotation as input
# for every 16S gene, identifies the closest 23S gene with the same seqid and strand values
# identifies unlined rRNA operons and those that cross the chromosome start/end boundary, removes these and writes their details to 'operons_unlinked_or_crossingboundary.tsv'
# writes coordinates of ITS regions to 'combined_ITS.gff'
# writes coordinates of 16S-ITS-23S regions full_operons.gff' for extraction later with bedtools
# writes 'operon_identifiers.txt' with the unique operon identifiers that are substituted into the fasta headers later
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/filter_make_gff_complete.R
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/filter_make_gff_incomplete.R

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
rm -f Outputs_Complete/combined_rrna.gff Outputs_Incomplete/combined_rrna.gff

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

# constructs database
mkdir -p Database/

# run dereplication
for i in 0.999 0.99 0.97 0.95 0.90
do

    # removes redundant sequences from Complete Assemblies dataset based on identity cutoffs
    # also extracts relevant cluster information
    vsearch \
        --cluster_smallmem Outputs_Complete/full_operons_sorted.fna \
        --id ${i} \
        --centroids Outputs_Complete/nrRep_${i}.fna \
        --uc Outputs_Complete/clusters_${i}.tsv \
        --threads ${i} \
        -log Outputs_Complete/vsearchlog_${i}.txt \
        --consout Outputs_Complete/nrCon_${i}.fna

    cat Outputs_Complete/clusters_${i}.tsv | \
        awk '$1=="S"' > Outputs_Complete/vsearch_centroids_${i}.tsv
    cat Outputs_Complete/clusters_${i}.tsv | \
        awk '$1=="C"' > Outputs_Complete/vsearch_clusters_${i}.tsv
    cat Outputs_Complete/clusters_${i}.tsv | \
        awk '$1=="H"' > Outputs_Complete/vsearch_hits_${i}.tsv

    # adds operon sequences from Incomplete Assemblies to the NR datasets generated from complete assemblies
    # then reruns clustering, and extracts relevant cluster information
    # only adds NR sequences from Incomplete genome assemblies if they are not already represented in the Complete dataset

    cat Outputs_Complete/nrRep_${i}.fna \
        Outputs_Incomplete/full_operons_sorted.fna > Outputs_Combined/full_operons_sorted_vsearchinput_${i}.fna

    vsearch \
        --cluster_smallmem Outputs_Combined/full_operons_sorted_vsearchinput_${i}.fna \
        --id ${i} \
        --centroids Outputs_Combined/nrRep_${i}.fna \
        --uc Outputs_Combined/clusters_${i}.tsv \
        --threads ${i} \
        -log Outputs_Combined/vsearchlog_${i}.txt \
        --consout Outputs_Combined/nrCon_${i}.fna \
        --usersort
        
    cat Outputs_Combined/clusters_${i}.tsv | \
        awk '$1=="S"' > Outputs_Combined/vsearch_centroids_${i}.tsv
    cat Outputs_Combined/clusters_${i}.tsv | \
        awk '$1=="C"' > Outputs_Combined/vsearch_clusters_${i}.tsv
    cat Outputs_Combined/clusters_${i}.tsv | \
        awk '$1=="H"' > Outputs_Combined/vsearch_hits_${i}.tsv

    rm -f Outputs_Combined/full_operons_sorted_vsearchinput_${i}.fna

    # edit fasta headers of nrCon sequence files
    sed -i 's/centroid=// ; s/;seqs=.*//' Outputs_Complete/nrCon_${i}.fna 
    sed -i 's/centroid=// ; s/;seqs=.*//' Outputs_Combined/nrCon_${i}.fna

    sortbyname.sh \
        in=Outputs_Combined/nrRep_${i}.fna \
        out=Database/nrRep_${i}.fna.gz

    sortbyname.sh \
        in=Outputs_Combined/nrCon_${i}.fna \
        out=Database/nrCon_${i}.fna.gz

    # remove database sequence files from output directories
    rm -f Outputs_Combined/nrRep_${i}.fna
    rm -f Outputs_Combined/nrCon_${i}.fna

done

# move Full database files to single directory
mv Outputs_Combined/full_operons_identifiers.fna Database/full.fna
gzip Database/full.fna

# get taxonomy for each rrn operon in databases
for i in 0.999 0.99 0.97 0.95 0.90
do

    Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_taxonomy_gtdb.R ${i}
    Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_taxonomy_ncbi.R ${i}

    cp Outputs_Combined/taxRep_gtdb_${i}.tsv Database/taxRep_gtdb_${i}.tsv
    cp Outputs_Combined/taxLCA_gtdb_${i}.tsv Database/taxLCA_gtdb_${i}.tsv
    cp Outputs_Combined/taxMaj_gtdb_${i}.tsv Database/taxMaj_gtdb_${i}.tsv

    cp Outputs_Combined/taxRep_ncbi_${i}.tsv Database/taxRep_ncbi_${i}.tsv
    cp Outputs_Combined/taxLCA_ncbi_${i}.tsv Database/taxLCA_ncbi_${i}.tsv
    cp Outputs_Combined/taxMaj_ncbi_${i}.tsv Database/taxMaj_ncbi_${i}.tsv

done

cp Outputs_Combined/taxFull_gtdb.tsv Database/taxFull_gtdb.tsv
cp Outputs_Combined/taxFull_ncbi.tsv Database/taxFull_ncbi.tsv

# calculate mean and median genome size for each taxon in the database (at all 7 taxonomic ranks)
# doing this using either all genomes or complete genomes only
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_genome_length_gtdb.R
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_genome_length_ncbi.R
cp stats_genomelength_gtdb.tsv stats_genomelength_ncbi.tsv Database/

# calculate mean and median rrn copy number for each taxon in the database (at all 7 taxonomic ranks)
# doing this using either all genomes or complete genomes only
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_copy_number_gtdb.R
Rscript ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_copy_number_ncbi.R
cp stats_copynumber_gtdb.tsv stats_copynumber_ncbi.tsv Database/

# move database to current working directory
mv Database/ ${VAR_SOURCE_DIRECTORY}/${VAR_OUTPUT_DIRECTORY}