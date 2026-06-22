#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

def main():
    if len(sys.argv) < 4:
        print("Usage: get_taxonomy.py <scheme: gtdb|ncbi> <seq_type: operon|16S|23S> <pid: identity>")
        sys.exit(1)
        
    scheme = sys.argv[1].lower()
    seq_type = sys.argv[2]
    pid = sys.argv[3]
    
    ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    # Load taxonomy.tsv (has header: accession, gtdb_taxonomy, ncbi_taxonomy)
    if not os.path.exists('taxonomy.tsv'):
        print("Error: taxonomy.tsv not found.", file=sys.stderr)
        sys.exit(1)
        
    tax_df = pd.read_csv('taxonomy.tsv', sep='\t')
    tax_df = tax_df[['accession', f'{scheme}_taxonomy']].copy()
    tax_df = tax_df.rename(columns={'accession': 'Assembly', f'{scheme}_taxonomy': 'Taxonomy'})
    tax_df['Taxonomy'] = tax_df['Taxonomy'].str.replace(';', '|')
    
    # Load operon mappings from master_rrna.gff
    if not os.path.exists('master_rrna.gff'):
        print("Error: master_rrna.gff not found.", file=sys.stderr)
        sys.exit(1)
        
    combined_operons = pd.read_csv('master_rrna.gff', sep='\t')
    combined_operons = combined_operons[['seqid', 'OperonID']].copy()
    combined_operons['Assembly'] = combined_operons['seqid'].apply(lambda x: x.split('__')[0])
    
    # Write taxFull file if it doesn't exist
    tax_full_name = f'taxFull_{seq_type}_{scheme}.tsv'
    if not os.path.exists(tax_full_name):
        tax_full = pd.merge(combined_operons, tax_df, on='Assembly', how='left')
        tax_full = tax_full[['OperonID', 'Taxonomy']].sort_values(by='OperonID')
        tax_full.to_csv(tax_full_name, sep='\t', header=False, index=False)
    
    # Load vsearch outputs
    centroids_file = f'vsearch_centroids_{seq_type}_{pid}.tsv'
    hits_file = f'vsearch_hits_{seq_type}_{pid}.tsv'
    
    if not os.path.exists(centroids_file):
        print(f"Centroids file {centroids_file} not found. Skipping.")
        sys.exit(0)
        
    try:
        centroids = pd.read_csv(centroids_file, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        centroids = pd.DataFrame(columns=range(10))
        
    try:
        hits = pd.read_csv(hits_file, sep='\t', header=None) if os.path.exists(hits_file) else pd.DataFrame(columns=range(10))
    except pd.errors.EmptyDataError:
        hits = pd.DataFrame(columns=range(10))
        
    if len(centroids) == 0:
        print(f"No centroids found in {centroids_file}.")
        sys.exit(0)
        
    # Map query sequence (column index 8) to Assembly via OperonID
    centroids_mapped = pd.merge(centroids[[1, 8]], combined_operons[['OperonID', 'Assembly']], left_on=8, right_on='OperonID')
    centroids_mapped = centroids_mapped.rename(columns={1: 'ClusterID'})
    
    if len(hits) > 0:
        hits_mapped = pd.merge(hits[[1, 8]], combined_operons[['OperonID', 'Assembly']], left_on=8, right_on='OperonID')
        hits_mapped = hits_mapped.rename(columns={1: 'ClusterID'})
        vsearch_df = pd.concat([centroids_mapped, hits_mapped], ignore_index=True)
    else:
        vsearch_df = centroids_mapped
        
    vsearch_df = vsearch_df.sort_values(by='ClusterID')
    
    # Split Taxonomy into individual ranks
    def split_taxonomy(tax_str):
        if pd.isna(tax_str):
            return ['NA'] * 7
        parts = tax_str.split('|')
        parts += ['NA'] * (7 - len(parts))
        return parts[:7]
        
    split_tax = tax_df['Taxonomy'].apply(split_taxonomy).tolist()
    tax_split_df = pd.DataFrame(split_tax, columns=ranks)
    tax_split_df['Assembly'] = tax_df['Assembly']
    tax_split_df['Taxonomy'] = tax_df['Taxonomy']
    
    vsearch_tax = pd.merge(vsearch_df, tax_split_df, on='Assembly', how='left')
    
    # --- 1. taxRep ---
    tax_rep_df = pd.merge(centroids_mapped[['ClusterID', 'OperonID', 'Assembly']], tax_df[['Assembly', 'Taxonomy']], on='Assembly', how='left')
    tax_rep_df = tax_rep_df[['OperonID', 'Taxonomy']].sort_values(by='OperonID')
    tax_rep_df.to_csv(f'taxRep_{seq_type}_{scheme}_{pid}.tsv', sep='\t', header=False, index=False)
    
    # --- 2. taxLCA ---
    lca_records = []
    for cluster_id, group in vsearch_tax.groupby('ClusterID'):
        lca_rank = 'Unknown'
        consensus_taxonomy = ['NA'] * 7
        for rank_idx, rank in enumerate(ranks):
            unique_vals = group[rank].dropna().unique()
            if len(unique_vals) == 1 and unique_vals[0] != 'NA' and unique_vals[0] != '':
                consensus_taxonomy[rank_idx] = unique_vals[0]
                lca_rank = rank
            else:
                break
        lca_records.append({
            'ClusterID': cluster_id,
            'LCA_Rank': lca_rank,
            **{rank: consensus_taxonomy[idx] for idx, rank in enumerate(ranks)}
        })
        
    lca_df = pd.DataFrame(lca_records)
    tax_lca = pd.merge(lca_df, centroids_mapped[['ClusterID', 'OperonID']], on='ClusterID')
    
    # LCA rank summary
    rank_summary = tax_lca['LCA_Rank'].value_counts().reset_index()
    rank_summary.columns = ['LCA_Rank', 'Count']
    rank_summary['Percentage'] = 100 * rank_summary['Count'] / rank_summary['Count'].sum()
    rank_summary.to_csv(f'taxLCA_{seq_type}_{scheme}_{pid}_ranksummary.tsv', sep='\t', index=False)
    
    def join_ranks(row):
        return '|'.join([row[r] for r in ranks])
        
    tax_lca['Taxonomy'] = tax_lca.apply(join_ranks, axis=1)
    tax_lca_out = tax_lca[['OperonID', 'Taxonomy']].sort_values(by='OperonID')
    tax_lca_out.to_csv(f'taxLCA_{seq_type}_{scheme}_{pid}.tsv', sep='\t', header=False, index=False)
    
    # --- 3. taxMaj ---
    maj_records = []
    for cluster_id, group in vsearch_tax.groupby('ClusterID'):
        maj_taxonomy = ['NA'] * 7
        for rank_idx, rank in enumerate(ranks):
            counts = group[rank].value_counts()
            if len(counts) > 0:
                max_count = counts.max()
                maj_vals = counts[counts == max_count].index.tolist()
                if len(maj_vals) == 1 and maj_vals[0] != 'NA' and maj_vals[0] != '':
                    maj_taxonomy[rank_idx] = maj_vals[0]
                else:
                    break
            else:
                break
        maj_records.append({
            'ClusterID': cluster_id,
            'Taxonomy': '|'.join(maj_taxonomy)
        })
        
    maj_df = pd.DataFrame(maj_records)
    tax_maj = pd.merge(maj_df, centroids_mapped[['ClusterID', 'OperonID']], on='ClusterID')
    
    def clean_maj_taxonomy(tax_str):
        parts = tax_str.split('|')
        last_valid = -1
        for idx, part in enumerate(parts):
            if part != 'NA' and part != '':
                last_valid = idx
            else:
                break
        if last_valid == -1:
            return 'NA'
        return '|'.join(parts[:last_valid+1])
        
    tax_maj['Taxonomy'] = tax_maj['Taxonomy'].apply(clean_maj_taxonomy)
    tax_maj_out = tax_maj[['OperonID', 'Taxonomy']].sort_values(by='OperonID')
    tax_maj_out.to_csv(f'taxMaj_{seq_type}_{scheme}_{pid}.tsv', sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()
