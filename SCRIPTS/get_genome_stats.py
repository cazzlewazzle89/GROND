#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

def expand_taxonomy(df, tax_col):
    rows = []
    for _, row in df.iterrows():
        tax_str = row[tax_col]
        if pd.isna(tax_str):
            continue
        # Remove strain suffix to match the R script
        tax_str = tax_str.split(';t__')[0]
        parts = tax_str.split(';')
        for i, part in enumerate(parts):
            if not part:
                continue
            rank_char = part.split('__')[0].upper()
            lineage = ';'.join(parts[:i+1])
            rows.append({
                'Assembly': row['Assembly'],
                'Value': row['Value'],
                'Taxonomy': part,
                'Rank': rank_char,
                'Lineage': lineage
            })
    return pd.DataFrame(rows)

def process_stats(df_expanded, value_name, complete_assemblies):
    if len(df_expanded) == 0:
        return pd.DataFrame()
        
    # 1. Combined (All) stats
    stats_all = df_expanded.groupby(['Taxonomy', 'Rank', 'Lineage'])['Value'].agg([
        ('Count', 'count'),
        ('Mean', 'mean'),
        ('Median', 'median'),
        ('SD', 'std')
    ]).reset_index()
    
    if value_name == 'CopyNumber':
        stats_all = stats_all.rename(columns={
            'Count': 'GenomeCount_All',
            'Mean': 'CopyNumber_Mean_All',
            'Median': 'CopyNumber_Median_All',
            'SD': 'CopyNumber_SD_All'
        })
    else:
        stats_all = stats_all.rename(columns={
            'Count': 'GenomeCount_All',
            'Mean': 'GenomeLength_Mean_All',
            'Median': 'GenomeLength_Median_All',
            'SD': 'GenomeLength_SD_All'
        })
        
    # 2. Complete stats
    df_complete = df_expanded[df_expanded['Assembly'].isin(complete_assemblies)]
    if len(df_complete) > 0:
        stats_complete = df_complete.groupby(['Taxonomy', 'Rank', 'Lineage'])['Value'].agg([
            ('Count', 'count'),
            ('Mean', 'mean'),
            ('Median', 'median'),
            ('SD', 'std')
        ]).reset_index()
    else:
        stats_complete = pd.DataFrame(columns=['Taxonomy', 'Rank', 'Lineage', 'Count', 'Mean', 'Median', 'SD'])
        
    if value_name == 'CopyNumber':
        stats_complete = stats_complete.rename(columns={
            'Count': 'GenomeCount_Complete',
            'Mean': 'CopyNumber_Mean_Complete',
            'Median': 'CopyNumber_Median_Complete',
            'SD': 'CopyNumber_SD_Complete'
        })
    else:
        stats_complete = stats_complete.rename(columns={
            'Count': 'GenomeCount_Complete',
            'Mean': 'GenomeLength_Mean_Complete',
            'Median': 'GenomeLength_Median_Complete',
            'SD': 'GenomeLength_SD_Complete'
        })
        
    # Merge them
    merged = pd.merge(stats_all, stats_complete, on=['Taxonomy', 'Rank', 'Lineage'], how='left')
    
    # Sort by rank
    rank_order = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5, 'S': 6}
    merged['Rank_Sort'] = merged['Rank'].map(rank_order)
    merged = merged.sort_values(by=['Rank_Sort', 'Taxonomy', 'Lineage']).drop(columns='Rank_Sort')
    
    return merged

def main():
    if len(sys.argv) < 3:
        print("Usage: get_genome_stats.py <scheme: gtdb|ncbi> <manifest_path>")
        sys.exit(1)
        
    scheme = sys.argv[1].lower()
    manifest_path = sys.argv[2]
    
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
    
    # Load taxonomy
    if not os.path.exists('taxonomy.tsv'):
        print("Error: taxonomy.tsv not found.", file=sys.stderr)
        sys.exit(1)
        
    tax_df = pd.read_csv('taxonomy.tsv', sep='\t')
    tax_df = tax_df[['accession', f'{scheme}_taxonomy']].copy()
    tax_df = tax_df.rename(columns={'accession': 'Assembly', f'{scheme}_taxonomy': 'Taxonomy'})
    
    # --- Part A: Genome Length Stats ---
    if not os.path.exists('Outputs/seq_length.tsv'):
        print("Error: Outputs/seq_length.tsv not found.", file=sys.stderr)
        sys.exit(1)
        
    seq_len_df = pd.read_csv('Outputs/seq_length.tsv', sep='\t', header=None, names=['Assembly__SeqID', 'Length'])
    seq_len_df['Assembly'] = seq_len_df['Assembly__SeqID'].apply(lambda x: x.split('__')[0])
    
    genome_len_df = seq_len_df.groupby('Assembly')['Length'].sum().reset_index()
    genome_len_df = pd.merge(genome_len_df, tax_df, on='Assembly', how='inner')
    genome_len_df = genome_len_df.rename(columns={'Length': 'Value'})
    
    expanded_len = expand_taxonomy(genome_len_df, 'Taxonomy')
    stats_len = process_stats(expanded_len, 'GenomeLength', complete_assemblies)
    stats_len.to_csv(f'stats_genomelength_{scheme}.tsv', sep='\t', index=False)
    
    # --- Part B: Copy Number Stats ---
    if not os.path.exists('Outputs/master_rrna.tsv'):
        print("Error: Outputs/master_rrna.tsv not found.", file=sys.stderr)
        sys.exit(1)
        
    master_rrna = pd.read_csv('Outputs/master_rrna.tsv', sep='\t')
    master_rrna['Assembly'] = master_rrna['seqid'].apply(lambda x: x.split('__')[0])
    
    copy_num_df = master_rrna.groupby('Assembly').size().reset_index(name='CopyNumber')
    copy_num_df = pd.merge(copy_num_df, tax_df, on='Assembly', how='inner')
    copy_num_df = copy_num_df.rename(columns={'CopyNumber': 'Value'})
    
    expanded_copy = expand_taxonomy(copy_num_df, 'Taxonomy')
    stats_copy = process_stats(expanded_copy, 'CopyNumber', complete_assemblies)
    stats_copy.to_csv(f'stats_copynumber_{scheme}.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()
