#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: filter_make_gff.py <manifest_path>")
        sys.exit(1)
        
    manifest_path = sys.argv[1]
    if not os.path.exists(manifest_path):
        print(f"Error: manifest file {manifest_path} not found.", file=sys.stderr)
        sys.exit(1)
        
    # Read manifest.tsv to identify complete assemblies
    manifest = pd.read_csv(manifest_path, sep='\t')
    manifest.columns = [c.strip() for c in manifest.columns]
    for col in manifest.columns:
        if manifest[col].dtype == object:
            manifest[col] = manifest[col].str.strip()
            
    complete_assemblies = set(manifest[manifest['genome_completeness'].str.lower() == 'complete']['genomeID'].tolist())

    # Read GFF. GFF is tab-separated, 9 columns.
    cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    try:
        df = pd.read_csv('combined_rrna.gff', sep='\t', comment='#', names=cols, header=None)
    except Exception as e:
        print(f"Error reading combined_rrna.gff: {e}", file=sys.stderr)
        sys.exit(1)

    # Sort to ensure relative ordering is correct
    df = df.sort_values(by=['seqid', 'start']).reset_index(drop=True)

    # Filter out partial genes
    df = df[~df['attributes'].str.contains('partial=true', case=False, na=False)].copy()

    # Extract product name
    def get_product(attr):
        match = re.search(r'product=([^;]+)', attr)
        if match:
            return match.group(1).strip()
        return ''

    df['product'] = df['attributes'].apply(get_product)

    # Filter for valid rRNA products
    valid_products = {
        '16S ribosomal RNA', '23S ribosomal RNA',
        '16S Ribosomal RNA', '23S Ribosomal RNA',
        'ribosomal RNA-16S', 'ribosomal RNA-23S'
    }
    df = df[df['product'].isin(valid_products)].copy()

    # Classify rRNA gene
    df['rrna_gene'] = df['product'].apply(lambda x: '16S' if '16S' in x else ('23S' if '23S' in x else 'Other'))
    df = df[df['rrna_gene'] != 'Other'].reset_index(drop=True)

    # Loop through and match 16S with the closest 23S on the same contig and strand
    rows_to_keep = []
    n = len(df)
    for i in range(n):
        row = df.iloc[i]
        if row['rrna_gene'] != '16S':
            continue

        if row['strand'] == '+':
            if i + 1 >= n:
                continue
            next_row = df.iloc[i + 1]
            if (next_row['rrna_gene'] == '23S' and 
                next_row['seqid'] == row['seqid'] and 
                next_row['strand'] == '+'):
                rows_to_keep.extend([i, i + 1])

        elif row['strand'] == '-':
            if i - 1 < 0:
                continue
            prev_row = df.iloc[i - 1]
            if (prev_row['rrna_gene'] == '23S' and 
                prev_row['seqid'] == row['seqid'] and 
                prev_row['strand'] == '-'):
                rows_to_keep.extend([i - 1, i])

    # De-duplicate indices
    rows_to_keep = sorted(list(set(rows_to_keep)))
    operons = df.iloc[rows_to_keep].copy().reset_index(drop=True)

    if len(operons) == 0:
        print("No valid operons identified.")
        sys.exit(0)

    # Assign unique IDs: {genomeID}__{original_seqid}__{num}
    # Since seqid is already prepended with genomeID__ (i.e. genomeID__original_seqid),
    # we just count within each seqid.
    operon_ids = []
    operon_counter = {}
    for k in range(len(operons) // 2):
        row = operons.iloc[2 * k]
        seqid = row['seqid']
        count = operon_counter.get(seqid, 0) + 1
        operon_counter[seqid] = count
        
        # Unique ID is seqid__count -> {genomeID}__{original_seqid}__{count}
        op_id = f"{seqid}__{count}"
        operon_ids.extend([op_id, op_id])

    operons['OperonID'] = operon_ids

    # Separate 16S and 23S to pivot wider
    df_16s = operons[operons['rrna_gene'] == '16S'].copy()
    df_23s = operons[operons['rrna_gene'] == '23S'].copy()

    df_16s = df_16s.rename(columns={'start': 'start_16S', 'end': 'end_16S'})
    df_23s = df_23s.rename(columns={'start': 'start_23S', 'end': 'end_23S'})

    merged = pd.merge(
        df_16s[['OperonID', 'seqid', 'strand', 'start_16S', 'end_16S']],
        df_23s[['OperonID', 'seqid', 'strand', 'start_23S', 'end_23S']],
        on=['OperonID', 'seqid', 'strand']
    )

    # Read sequence lengths
    try:
        seq_len_df = pd.read_csv('seq_length.tsv', sep='\t', header=None, names=['seqid', 'Sequence_Length'])
    except Exception as e:
        print(f"Error reading seq_length.tsv: {e}", file=sys.stderr)
        sys.exit(1)

    merged = pd.merge(merged, seq_len_df, on='seqid', how='left')

    # Calculate ITS Length
    merged['ITS_Length'] = np.where(
        merged['strand'] == '+',
        merged['start_23S'] - merged['end_16S'],
        merged['start_16S'] - merged['end_23S']
    )

    # Determine Linkage
    merged['Linkage'] = np.where(merged['ITS_Length'] > 1500, 'Unlinked', 'Linked')

    # Calculate boundaries
    coords = merged[['start_16S', 'end_16S', 'start_23S', 'end_23S']]
    merged['Operon_Start'] = coords.min(axis=1)
    merged['Operon_End'] = coords.max(axis=1)

    merged['CrossBoundary'] = (merged['Operon_Start'] == 1) | (merged['Operon_End'] >= merged['Sequence_Length'])

    # Write unlinked/boundary-crossing operons to file
    unlinked_boundary = merged[(merged['CrossBoundary'] == True) | (merged['Linkage'] == 'Unlinked')]
    unlinked_boundary.to_csv('operons_unlinked_or_crossingboundary.tsv', sep='\t', index=False)

    # Filter for good operons
    good_operons = merged[(merged['CrossBoundary'] == False) & (merged['Linkage'] == 'Linked')].copy()

    good_operons['source'] = 'Dummy'
    good_operons['type'] = 'rRNA'
    good_operons['score'] = '.'
    good_operons['phase'] = '.'

    good_operons['start_ITS'] = np.where(
        good_operons['strand'] == '+',
        good_operons['end_16S'],
        good_operons['end_23S']
    )
    good_operons['end_ITS'] = np.where(
        good_operons['strand'] == '+',
        good_operons['start_23S'],
        good_operons['start_16S']
    )

    # Write master_rrna.gff
    master_cols = [
        'seqid', 'OperonID', 'source', 'type', 'start_ITS', 'end_ITS',
        'start_16S', 'end_16S', 'start_23S', 'end_23S', 'score', 'strand', 'phase'
    ]
    good_operons[master_cols].to_csv('master_rrna.gff', sep='\t', index=False)

    # Write full_ITS.gff
    its_gff = pd.DataFrame({
        'seqid_operon': good_operons['seqid'] + ' ' + good_operons['OperonID'],
        'source': good_operons['source'],
        'type': good_operons['type'],
        'start': good_operons['start_ITS'],
        'end': good_operons['end_ITS'],
        'score': good_operons['score'],
        'strand': good_operons['strand'],
        'phase': good_operons['phase']
    })
    its_gff.to_csv('full_ITS.gff', sep='\t', header=False, index=False)

    # Helper to check if a seqid belongs to a complete genome
    def is_complete_seqid(seqid):
        genome_id = seqid.split('__')[0]
        return genome_id in complete_assemblies

    # Split into complete and incomplete subsets
    complete_mask = good_operons['seqid'].apply(is_complete_seqid)
    operons_complete = good_operons[complete_mask].copy()
    operons_incomplete = good_operons[~complete_mask].copy()

    # Write separate coordinate files for bedtools getfasta
    def write_subset_gffs(df, prefix):
        if len(df) == 0:
            # Create empty files
            open(f'operon_{prefix}.gff', 'w').close()
            open(f'16S_{prefix}.gff', 'w').close()
            open(f'23S_{prefix}.gff', 'w').close()
            open(f'operon_identifiers_{prefix}.txt', 'w').close()
            return

        # Operons GFF
        pd.DataFrame({
            'seqid': df['seqid'],
            'source': df['source'],
            'type': df['type'],
            'start': df['Operon_Start'],
            'end': df['Operon_End'],
            'score': df['score'],
            'strand': df['strand'],
            'phase': df['phase']
        }).to_csv(f'operon_{prefix}.gff', sep='\t', header=False, index=False)

        # 16S GFF
        pd.DataFrame({
            'seqid': df['seqid'],
            'source': df['source'],
            'type': df['type'],
            'start': df[['start_16S', 'end_16S']].min(axis=1),
            'end': df[['start_16S', 'end_16S']].max(axis=1),
            'score': df['score'],
            'strand': df['strand'],
            'phase': df['phase']
        }).to_csv(f'16S_{prefix}.gff', sep='\t', header=False, index=False)

        # 23S GFF
        pd.DataFrame({
            'seqid': df['seqid'],
            'source': df['source'],
            'type': df['type'],
            'start': df[['start_23S', 'end_23S']].min(axis=1),
            'end': df[['start_23S', 'end_23S']].max(axis=1),
            'score': df['score'],
            'strand': df['strand'],
            'phase': df['phase']
        }).to_csv(f'23S_{prefix}.gff', sep='\t', header=False, index=False)

        # Identifiers mapping file (col 1: seqid, col 2: OperonID)
        df[['seqid', 'OperonID']].to_csv(f'operon_identifiers_{prefix}.txt', sep='\t', header=False, index=False)

    write_subset_gffs(operons_complete, 'complete')
    write_subset_gffs(operons_incomplete, 'incomplete')

if __name__ == '__main__':
    main()
