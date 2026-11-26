#!/usr/bin/env bash

# Check a GROND manifest for genomes whose FASTA assembly exists but whose barrnap GFF does not.
# For each such genome, run barrnap to create the GFF in a local directory, then write an updated
# manifest with path_to_gff pointing at the new file. The input manifest is never modified.
#
# Genomes with a missing FASTA are left unchanged in the updated manifest (grond.sh will report them).
# Already-generated GFFs in the output directory are reused, so an interrupted run can be resumed.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 -i <manifest.tsv> [-o <updated_manifest.tsv>] [-d <gff_dir>] [-t <threads_per_job>] [-j <parallel_jobs>] [-T <tmp_dir>]

  -i  Input manifest (genomeID, path_to_fna, path_to_gff, gtdb_taxonomy, ncbi_taxonomy, genome_completeness)
  -o  Updated manifest to write          (default: <input_basename>_updated.tsv in the current directory)
  -d  Directory for generated GFF files  (default: ./barrnap_gff)
  -t  Threads per barrnap job            (default: 1)
  -j  Number of barrnap jobs in parallel (default: 1)
  -T  Directory for temporary files      (default: \${TMPDIR:-/tmp})
EOF
    exit 1
}

MANIFEST_FILE=""
OUT_MANIFEST=""
GFF_DIR="${PWD}/barrnap_gff"
THREADS=1
JOBS=1
TMP_BASE="${TMPDIR:-/tmp}"

while getopts ":i:o:d:t:j:T:h" opt; do
    case $opt in
        i) MANIFEST_FILE="$OPTARG" ;;
        o) OUT_MANIFEST="$OPTARG" ;;
        d) GFF_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        T) TMP_BASE="$OPTARG" ;;
        *) usage ;;
    esac
done

[ -z "$MANIFEST_FILE" ] && usage
if [ ! -f "$MANIFEST_FILE" ]; then
    echo "ERROR: manifest not found: $MANIFEST_FILE" >&2
    exit 1
fi
if [ -z "$OUT_MANIFEST" ]; then
    base=$(basename "$MANIFEST_FILE")
    OUT_MANIFEST="${PWD}/${base%.*}_updated.tsv"
fi

# Resolve to absolute paths
abspath() { [[ "$1" = /* ]] && echo "$1" || echo "${PWD}/$1"; }
MANIFEST_FILE=$(abspath "$MANIFEST_FILE")
OUT_MANIFEST=$(abspath "$OUT_MANIFEST")
GFF_DIR=$(abspath "$GFF_DIR")

if [ "$MANIFEST_FILE" = "$OUT_MANIFEST" ]; then
    echo "ERROR: output manifest must differ from the input manifest." >&2
    exit 1
fi

if ! command -v barrnap >/dev/null 2>&1; then
    echo "ERROR: barrnap not found in \$PATH (activate the grond conda environment)." >&2
    exit 1
fi

VAR_WORKING_DIRECTORY=${PWD}
mkdir -p "$GFF_DIR"
WORK_DIR=$(mktemp -d "${TMP_BASE}/fill_missing_gff.XXXXXX")
trap 'rm -rf "$WORK_DIR"' EXIT

TODO_FILE="${WORK_DIR}/todo.tsv"
FAILED_FILE="${GFF_DIR}/failed_genomes.tsv"
: > "$TODO_FILE"
echo -e "genomeID\tpath_to_fna\treason" > "$FAILED_FILE"

resolve() {
    # Resolve manifest paths relative to the invocation directory, as grond.sh does
    local p="$1"
    if [ -n "$p" ] && [[ ! "$p" = /* ]]; then
        p="${VAR_WORKING_DIRECTORY}/${p}"
    fi
    echo "$p"
}

clean() { echo "$1" | tr -d '\r' | xargs; }

# --- Pass 1: find genomes with a FASTA but no GFF ---
n_total=0; n_ok=0; n_missing_fna=0; n_todo=0
while IFS=$'\t' read -r genomeID path_to_fna path_to_gff gtdb_taxonomy rest || [ -n "$genomeID" ]; do
    genomeID=$(clean "$genomeID")
    if [ "$genomeID" = "genomeID" ] || [ -z "$genomeID" ] || [[ "$genomeID" =~ ^# ]]; then
        continue
    fi
    n_total=$((n_total + 1))
    fna=$(resolve "$(clean "$path_to_fna")")
    gff=$(resolve "$(clean "$path_to_gff")")

    if [ ! -f "$fna" ]; then
        n_missing_fna=$((n_missing_fna + 1))
        continue
    fi
    if [ -n "$(clean "$path_to_gff")" ] && [ -f "$gff" ]; then
        n_ok=$((n_ok + 1))
        continue
    fi

    # Choose barrnap kingdom from the GTDB domain (default bacteria)
    kingdom=bac
    [[ "$gtdb_taxonomy" == *d__Archaea* ]] && kingdom=arc

    echo -e "${genomeID}\t${fna}\t${kingdom}" >> "$TODO_FILE"
    n_todo=$((n_todo + 1))
done < "$MANIFEST_FILE"

echo "Manifest checked: ${n_total} genomes; ${n_ok} with FASTA and GFF; ${n_missing_fna} with missing FASTA; ${n_todo} need a GFF."

# --- Run barrnap on genomes that need a GFF ---
run_barrnap() {
    local genomeID="$1" fna="$2" kingdom="$3"
    local out="${GFF_DIR}/${genomeID}_rrna.gff"
    [ -s "$out" ] && return 0   # reuse from a previous run

    local input="$fna"
    if [[ "$fna" == *.gz ]]; then
        input="${WORK_DIR}/${genomeID}.fna"
        if ! gunzip -c "$fna" > "$input"; then
            echo -e "${genomeID}\t${fna}\tdecompression failed" >> "$FAILED_FILE"
            rm -f "$input"; return 0
        fi
    fi

    # Write to a temporary file and move into place, so partial outputs are never mistaken for complete ones
    if barrnap --threads "$THREADS" --kingdom "$kingdom" --quiet "$input" > "${out}.tmp" 2> "${WORK_DIR}/${genomeID}.log"; then
        mv "${out}.tmp" "$out"
    else
        echo -e "${genomeID}\t${fna}\tbarrnap failed: $(tail -n 1 "${WORK_DIR}/${genomeID}.log")" >> "$FAILED_FILE"
        rm -f "${out}.tmp"
    fi
    [ "$input" != "$fna" ] && rm -f "$input"
    return 0
}
export -f run_barrnap
export GFF_DIR WORK_DIR FAILED_FILE THREADS

if [ "$n_todo" -gt 0 ]; then
    echo "Running barrnap (kingdom from GTDB domain) on ${n_todo} genomes with ${JOBS} parallel job(s) x ${THREADS} thread(s)..."
    tr '\t\n' '\0\0' < "$TODO_FILE" \
        | xargs -0 -n 3 -P "$JOBS" bash -c 'run_barrnap "$@"' _
fi

# --- Pass 2: write updated manifest ---
# Only genomes whose GFF now exists get a new path; everything else is copied through unchanged
DONE_FILE="${WORK_DIR}/done.txt"
: > "$DONE_FILE"
while IFS=$'\t' read -r genomeID _ _; do
    [ -s "${GFF_DIR}/${genomeID}_rrna.gff" ] && echo "$genomeID" >> "$DONE_FILE"
done < "$TODO_FILE"

awk -v gffdir="$GFF_DIR" 'BEGIN{FS=OFS="\t"}
    NR==FNR {done[$1]=1; next}
    {
        id=$1; gsub(/\r|^[ \t]+|[ \t]+$/, "", id)
        if (id in done) $3 = gffdir "/" id "_rrna.gff"
        print
    }' "$DONE_FILE" "$MANIFEST_FILE" > "$OUT_MANIFEST"
n_filled=$(wc -l < "$DONE_FILE" | tr -d ' ')

n_failed=$(($(wc -l < "$FAILED_FILE") - 1))
echo "Done: ${n_filled} GFFs created/reused in ${GFF_DIR}; ${n_failed} failed."
echo "Updated manifest written to ${OUT_MANIFEST}"
if [ "$n_failed" -gt 0 ]; then
    echo "WARNING: failed genomes listed in ${FAILED_FILE} (their manifest rows are unchanged)." >&2
else
    rm -f "$FAILED_FILE"
fi
