#!/usr/bin/env python3
"""Map every sequence in the full (pre-dereplication) dataset to its final
representative sequence at each clustering threshold.

Clustering in grond.sh is two-stage:
  1. complete-genome sequences are clustered alone  -> clusters_complete_{type}_{pid}.tsv
  2. stage-1 centroids + incomplete-genome sequences -> clusters_{type}_{pid}.tsv
A complete-genome sequence absorbed in stage 1 therefore never appears in the
stage-2 .uc file, so its final representative is resolved through its stage-1
centroid: final[x] = stage2[stage1[x]].

Output columns: OperonID, Source (complete|incomplete), then one column per
threshold (Rep_<pid>) holding the OperonID of the representative sequence in
<type>_nrRep_<pid>.fna.gz.
"""
import sys
import os
import pandas as pd


def read_uc(path):
    """Return {query_label: centroid_label} from the S and H records of a .uc file."""
    mapping = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return mapping
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts[0] == 'S':
                mapping[parts[8]] = parts[8]
            elif parts[0] == 'H':
                mapping[parts[8]] = parts[9]
    return mapping


def read_headers(path):
    ids = []
    if not os.path.exists(path):
        return ids
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def main():
    if len(sys.argv) < 4:
        print("Usage: make_cluster_map.py <seq_type: operon|16S|23S> <output_tsv> <pid> [<pid> ...]")
        sys.exit(1)

    seq_type = sys.argv[1]
    out_file = sys.argv[2]
    pids = sys.argv[3:]

    complete_ids = read_headers(f'Outputs/{seq_type}_complete_sorted.fna')
    incomplete_ids = read_headers(f'Outputs/{seq_type}_incomplete_sorted.fna')

    df = pd.DataFrame({
        'OperonID': complete_ids + incomplete_ids,
        'Source': ['complete'] * len(complete_ids) + ['incomplete'] * len(incomplete_ids),
    })

    n_missing_total = 0
    for pid in pids:
        stage1 = read_uc(f'Outputs/clusters_complete_{seq_type}_{pid}.tsv')
        stage2 = read_uc(f'Outputs/clusters_{seq_type}_{pid}.tsv')

        def resolve(seq_id):
            return stage2.get(stage1.get(seq_id, seq_id), pd.NA)

        df[f'Rep_{pid}'] = df['OperonID'].map(resolve)
        n_missing = int(df[f'Rep_{pid}'].isna().sum())
        n_reps = df[f'Rep_{pid}'].nunique()
        print(f"[{seq_type}/{pid}] {len(df)} sequences -> {n_reps} representatives"
              + (f"; WARNING: {n_missing} sequences unresolved" if n_missing else ""))
        n_missing_total += n_missing

    df = df.sort_values(by='OperonID')
    df.to_csv(out_file, sep='\t', index=False, na_rep='NA')

    if n_missing_total:
        print(f"WARNING [{seq_type}]: {n_missing_total} sequence/threshold combinations "
              f"could not be mapped to a representative.", file=sys.stderr)


if __name__ == '__main__':
    main()
