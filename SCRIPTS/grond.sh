#!/usr/bin/env bash

# Exit immediately if any command exits with a non-zero status
set -e

if [ -z "$1" ]; then
    echo "Usage: $0 <manifest.tsv>"
    exit 1
fi

# Resolve manifest path to absolute path
if [[ "$1" = /* ]]; then
    MANIFEST_FILE="$1"
else
    MANIFEST_FILE="${PWD}/$1"
fi

VAR_THREADS_SEQLENGTH=24
VAR_THREADS_VSEARCH=24
VAR_OUTPUT_DIRECTORY=R232
VAR_WORKING_DIRECTORY=${PWD}
VAR_TEMP_DIRECTORY=${PWD}/TEMP

VAR_SOURCE_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Create temp directory and output directory structure
mkdir -p ${VAR_TEMP_DIRECTORY}/Outputs
mkdir -p ${VAR_TEMP_DIRECTORY}/Database
cd ${VAR_TEMP_DIRECTORY}

echo "Parsing manifest file..."
# Prepare combined files
cat /dev/null > Outputs/combined.fna
cat /dev/null > Outputs/combined_rrna.gff
echo -e "accession\tgtdb_taxonomy\tncbi_taxonomy" > taxonomy.tsv

# Read manifest.tsv line by line, skipping header
while IFS=$'\t' read -r genomeID path_to_fna path_to_gff gtdb_taxonomy ncbi_taxonomy genome_completeness || [ -n "$genomeID" ]; do
    # Strip carriage returns and spaces
    genomeID=$(echo "$genomeID" | tr -d '\r' | xargs)
    path_to_fna=$(echo "$path_to_fna" | tr -d '\r' | xargs)
    path_to_gff=$(echo "$path_to_gff" | tr -d '\r' | xargs)
    gtdb_taxonomy=$(echo "$gtdb_taxonomy" | tr -d '\r' | xargs)
    ncbi_taxonomy=$(echo "$ncbi_taxonomy" | tr -d '\r' | xargs)
    
    # Skip header row or empty lines
    if [ "$genomeID" = "genomeID" ] || [ -z "$genomeID" ] || [[ "$genomeID" =~ ^# ]]; then
        continue
    fi
    
    # Resolve relative paths against the invocation directory, not the software directory
    if [[ ! "$path_to_fna" = /* ]]; then
        path_to_fna="${VAR_WORKING_DIRECTORY}/${path_to_fna}"
    fi
    if [[ ! "$path_to_gff" = /* ]]; then
        path_to_gff="${VAR_WORKING_DIRECTORY}/${path_to_gff}"
    fi
    
    missing=0
    if [ ! -f "$path_to_fna" ]; then
        echo "WARNING: Skipping ${genomeID} — FASTA file not found: $path_to_fna"
        missing=1
    fi
    if [ ! -f "$path_to_gff" ]; then
        echo "WARNING: Skipping ${genomeID} — GFF file not found: $path_to_gff"
        missing=1
    fi
    [ "$missing" -eq 1 ] && continue
    
    # 1. Edit fasta headers to contain the genome ID
    if [[ "$path_to_fna" == *.gz ]]; then
        gunzip -c "$path_to_fna"
    else
        cat "$path_to_fna"
    fi | sed "s/^>/>${genomeID}__/ ; s/ .*//" >> Outputs/combined.fna
    
    # 2. Edit GFF seqid column to contain the genome ID
    if [[ "$path_to_gff" == *.gz ]]; then
        gunzip -c "$path_to_gff"
    else
        cat "$path_to_gff"
    fi | awk -v prefix="${genomeID}__" 'BEGIN{FS=OFS="\t"} !/^#/ {$1=prefix $1} {print}' >> Outputs/combined_rrna.gff
    
    # 3. Add to taxonomy.tsv
    echo -e "${genomeID}\t${gtdb_taxonomy}\t${ncbi_taxonomy}" >> taxonomy.tsv

done < "$MANIFEST_FILE"

echo "Calculating sequence lengths with seqkit..."
seqkit fx2tab Outputs/combined.fna -n -l -j ${VAR_THREADS_SEQLENGTH} > Outputs/seq_length.tsv

echo "Running GFF filtering and coordinate extraction..."
python3 ${VAR_SOURCE_DIRECTORY}/SCRIPTS/filter_make_gff.py "$MANIFEST_FILE"

echo "Extracting and sorting fasta sequences..."
for type in operon 16S 23S; do
    for comp in complete incomplete; do
        if [ -s Outputs/${type}_${comp}.gff ]; then
            bedtools getfasta \
                -fi Outputs/combined.fna \
                -bed Outputs/${type}_${comp}.gff \
                -fo Outputs/${type}_${comp}.fna \
                -name \
                -s
            python3 ${VAR_SOURCE_DIRECTORY}/SCRIPTS/replace_headers.py \
                Outputs/${type}_${comp}.fna \
                Outputs/operon_identifiers_${comp}.txt \
                Outputs/${type}_${comp}_identifiers.fna
            sortbyname.sh \
                in=Outputs/${type}_${comp}_identifiers.fna \
                out=Outputs/${type}_${comp}_sorted.fna \
                length \
                descending
        else
            touch Outputs/${type}_${comp}_sorted.fna
        fi
    done
done

# Prioritized clustering function
run_clustering() {
    local type=$1
    local pid=$2
    
    echo "Clustering $type at identity threshold $pid..."
    
    # 1. Cluster complete genomes
    vsearch \
        --cluster_smallmem Outputs/${type}_complete_sorted.fna \
        --id ${pid} \
        --centroids Outputs/vsearch_centroids_raw_${type}_${pid}.fna \
        --uc Outputs/clusters_complete_${type}_${pid}.tsv \
        --threads ${VAR_THREADS_VSEARCH} \
        -log Outputs/vsearchlog_complete_${type}_${pid}.txt \
        --consout Outputs/vsearch_consensus_raw_${type}_${pid}.fna
        
    # 2. Combine complete centroids and incomplete genomes
    cat Outputs/vsearch_centroids_raw_${type}_${pid}.fna \
        Outputs/${type}_incomplete_sorted.fna > Outputs/vsearch_input_${type}_${pid}.fna
        
    # 3. Cluster combined dataset with --usersort (maintains priority)
    vsearch \
        --cluster_smallmem Outputs/vsearch_input_${type}_${pid}.fna \
        --id ${pid} \
        --centroids Outputs/vsearch_centroids_${type}_${pid}.fna \
        --uc Outputs/clusters_${type}_${pid}.tsv \
        --threads ${VAR_THREADS_VSEARCH} \
        -log Outputs/vsearchlog_${type}_${pid}.txt \
        --consout Outputs/vsearch_consensus_${type}_${pid}.fna \
        --usersort
        
    # 4. Extract centroids, hits, clusters
    cat Outputs/clusters_${type}_${pid}.tsv | awk '$1=="S"' > Outputs/vsearch_centroids_${type}_${pid}.tsv
    cat Outputs/clusters_${type}_${pid}.tsv | awk '$1=="C"' > Outputs/vsearch_clusters_${type}_${pid}.tsv
    cat Outputs/clusters_${type}_${pid}.tsv | awk '$1=="H"' > Outputs/vsearch_hits_${type}_${pid}.tsv
    
    # 5. Clean up consensus headers if consensus file is generated
    if [ -f Outputs/vsearch_consensus_${type}_${pid}.fna ]; then
        sed 's/centroid=// ; s/;seqs=.*//' Outputs/vsearch_consensus_${type}_${pid}.fna > Outputs/vsearch_consensus_${type}_${pid}.tmp && \
        mv Outputs/vsearch_consensus_${type}_${pid}.tmp Outputs/vsearch_consensus_${type}_${pid}.fna
    fi
    
    # 6. Sort and zip representative databases
    if [ "$type" = "operon" ]; then
        sortbyname.sh \
            in=Outputs/vsearch_centroids_${type}_${pid}.fna \
            out=Database/nrRep_${pid}.fna.gz
        if [ -f Outputs/vsearch_consensus_${type}_${pid}.fna ]; then
            sortbyname.sh \
                in=Outputs/vsearch_consensus_${type}_${pid}.fna \
                out=Database/nrCon_${pid}.fna.gz
        fi
    else
        sortbyname.sh \
            in=Outputs/vsearch_centroids_${type}_${pid}.fna \
            out=Database/${type}_nrRep_${pid}.fna.gz
        if [ -f Outputs/vsearch_consensus_${type}_${pid}.fna ]; then
            sortbyname.sh \
                in=Outputs/vsearch_consensus_${type}_${pid}.fna \
                out=Database/${type}_nrCon_${pid}.fna.gz
        fi
    fi
    
    # 7. Run taxonomy assignment python scripts
    python3 ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_taxonomy.py gtdb ${type} ${pid}
    python3 ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_taxonomy.py ncbi ${type} ${pid}
    
    # Copy taxonomy outputs to Database/
    if [ "$type" = "operon" ]; then
        cp Outputs/taxRep_operon_gtdb_${pid}.tsv Database/taxRep_gtdb_${pid}.tsv
        cp Outputs/taxLCA_operon_gtdb_${pid}.tsv Database/taxLCA_gtdb_${pid}.tsv
        cp Outputs/taxMaj_operon_gtdb_${pid}.tsv Database/taxMaj_gtdb_${pid}.tsv
        
        cp Outputs/taxRep_operon_ncbi_${pid}.tsv Database/taxRep_ncbi_${pid}.tsv
        cp Outputs/taxLCA_operon_ncbi_${pid}.tsv Database/taxLCA_ncbi_${pid}.tsv
        cp Outputs/taxMaj_operon_ncbi_${pid}.tsv Database/taxMaj_ncbi_${pid}.tsv
    else
        cp Outputs/taxRep_${type}_gtdb_${pid}.tsv Database/${type}_taxRep_gtdb_${pid}.tsv
        cp Outputs/taxLCA_${type}_gtdb_${pid}.tsv Database/${type}_taxLCA_gtdb_${pid}.tsv
        cp Outputs/taxMaj_${type}_gtdb_${pid}.tsv Database/${type}_taxMaj_gtdb_${pid}.tsv
        
        cp Outputs/taxRep_${type}_ncbi_${pid}.tsv Database/${type}_taxRep_ncbi_${pid}.tsv
        cp Outputs/taxLCA_${type}_ncbi_${pid}.tsv Database/${type}_taxLCA_ncbi_${pid}.tsv
        cp Outputs/taxMaj_${type}_ncbi_${pid}.tsv Database/${type}_taxMaj_ncbi_${pid}.tsv
    fi
    
    # Clean up intermediate files
    rm -f Outputs/vsearch_centroids_raw_${type}_${pid}.fna
    rm -f Outputs/vsearch_consensus_raw_${type}_${pid}.fna
    rm -f Outputs/vsearch_input_${type}_${pid}.fna
    rm -f Outputs/vsearch_centroids_${type}_${pid}.fna
    rm -f Outputs/vsearch_consensus_${type}_${pid}.fna
}

# Run operon clustering
for pid in 1.0 0.999 0.99 0.97 0.95 0.90; do
    run_clustering operon ${pid}
done

# Run 16S clustering
for pid in 1.0 0.99 0.98 0.97; do
    run_clustering 16S ${pid}
done

# Run 23S clustering
for pid in 1.0 0.99 0.98 0.97; do
    run_clustering 23S ${pid}
done

echo "Copying taxFull files to Database..."
cp Outputs/taxFull_operon_gtdb.tsv Database/taxFull_gtdb.tsv
cp Outputs/taxFull_operon_ncbi.tsv Database/taxFull_ncbi.tsv
cp Outputs/taxFull_16S_gtdb.tsv Database/16S_taxFull_gtdb.tsv
cp Outputs/taxFull_16S_ncbi.tsv Database/16S_taxFull_ncbi.tsv
cp Outputs/taxFull_23S_gtdb.tsv Database/23S_taxFull_gtdb.tsv
cp Outputs/taxFull_23S_ncbi.tsv Database/23S_taxFull_ncbi.tsv

echo "Calculating genome stats..."
python3 ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_genome_stats.py gtdb "$MANIFEST_FILE"
python3 ${VAR_SOURCE_DIRECTORY}/SCRIPTS/get_genome_stats.py ncbi "$MANIFEST_FILE"

cp stats_genomelength_gtdb.tsv Database/stats_genomelength_gtdb.tsv
cp stats_genomelength_ncbi.tsv Database/stats_genomelength_ncbi.tsv
cp stats_copynumber_gtdb.tsv Database/stats_copynumber_gtdb.tsv
cp stats_copynumber_ncbi.tsv Database/stats_copynumber_ncbi.tsv

echo "Copying master GFF to Database..."
cp Outputs/master_rrna.tsv Database/master_rrna.tsv

# Copy database to output directory (relative to invocation directory, outside TEMP)
mkdir -p ${VAR_WORKING_DIRECTORY}/${VAR_OUTPUT_DIRECTORY}
cp -r Database/* ${VAR_WORKING_DIRECTORY}/${VAR_OUTPUT_DIRECTORY}/

# Cleanup temp directory
rm -rf ${VAR_TEMP_DIRECTORY}

echo "GROND build completed successfully."