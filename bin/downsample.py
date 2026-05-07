#!/usr/bin/env python3
"""
Downsample single-cell RNA-seq counts using scanpy's multinomial sampling.
"""

import argparse
import os
from weakref import ref
import numpy as np
import scanpy as sc
import anndata as ad
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Downsample scRNA-seq counts to target depth"
    )
    parser.add_argument(
        '--input', type=str, default="/space/grp/rschwartz/rschwartz/downsampling-experiment/data/homo_sapiens/Human_Multiple_Cortical_Areas_SMART-seq.h5ad",
        help='Input h5ad file with raw counts'
    )
    parser.add_argument(
        '--output', type=str, default="/space/grp/rschwartz/rschwartz/downsampling-experiment/data/homo_sapiens/Human_Multiple_Cortical_Areas_SMART-seq_downsampled_10000.h5ad",
        help='Output h5ad file'
    )
    parser.add_argument(
        '--target_depth', type=float, default=10000,
        help='Target number of reads per cell'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '--counts_layer', type=str, default=None,
        help='Layer containing raw counts (default: use .X)'
    )
    args, _ = parser.parse_known_args()
    return args


def main():
    args = parse_arguments()

    print(f"Loading {args.input}...")
    adata = ad.read_h5ad(args.input)
    print(f"Dataset shape: {adata.shape}")

    # If counts are in a layer, move them to .X
    if args.counts_layer:
        print(f"Using counts from layer '{args.counts_layer}'")
        adata.X = adata.layers[args.counts_layer]

    # Store original counts
    adata.layers['original_counts'] = adata.X.copy()
    original_counts = np.array(adata.X.sum(axis=1)).flatten()
    print(f"Original depth - Mean: {original_counts.mean():.0f}, Median: {np.median(original_counts):.0f}")

    # Downsample using scanpy
    print(f"\nDownsampling to {args.target_depth} reads per cell (seed={args.seed})...")
    # Per-cell multinomial downsampling to target depth


# doing this wrong
# 
    np.random.seed(args.seed)
    X = adata.X
    if not isinstance(X, np.ndarray):
        X = X.toarray()

    target = int(args.target_depth)
    new_X = np.zeros_like(X, dtype=np.float64)
    for i in range(X.shape[0]): # for each cell
        cell_counts = X[i].astype(np.float64)  # vector of gene counts for this cell
        total = cell_counts.sum() # total counts for this cell

        if total <= target:
            new_X[i] = cell_counts
        elif total > 0:
            p = cell_counts / total  # vector gene probabilities for multinomial sampling
            new_X[i] = np.random.multinomial(target, p)

    adata.X = new_X

    # Record stats
    new_counts = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['downsampled_total_counts'] = new_counts
    adata.obs['original_total_counts'] = original_counts
    print(f"Downsampled depth - Mean: {new_counts.mean():.0f}, Median: {np.median(new_counts):.0f}")

    # Compute gene detection stats
    # Genes detected per cell (count > 0)
    genes_per_cell = (new_X > 0).sum(axis=1)
    avg_genes_per_cell = genes_per_cell.mean()
    median_genes_per_cell = np.median(genes_per_cell)

    # Total genes detected in dataset (at least one cell has count > 0)
    total_genes_detected = (new_X.sum(axis=0) > 0).sum()
    total_genes = new_X.shape[1]

    # Same stats for original data
    orig_genes_per_cell = (X > 0).sum(axis=1)
    orig_avg_genes_per_cell = orig_genes_per_cell.mean()
    orig_median_genes_per_cell = np.median(orig_genes_per_cell)
    orig_total_genes_detected = (X.sum(axis=0) > 0).sum()

    # Save gene detection stats to TSV
    import pandas as pd
    stats_df = pd.DataFrame([{
        'target_depth': args.target_depth,
        'seed': args.seed,
        'n_cells': new_X.shape[0],
        'total_genes': total_genes,
        'original_avg_counts_per_cell': original_counts.mean(),
        'downsampled_avg_counts_per_cell': new_counts.mean(),
        'original_avg_genes_per_cell': orig_avg_genes_per_cell,
        'downsampled_avg_genes_per_cell': avg_genes_per_cell,
        'original_median_genes_per_cell': orig_median_genes_per_cell,
        'downsampled_median_genes_per_cell': median_genes_per_cell,
        'original_total_genes_detected': orig_total_genes_detected,
        'downsampled_total_genes_detected': total_genes_detected
    }])
    stats_output = args.output.replace('.h5ad', '_stats.tsv')
    stats_df.to_csv(stats_output, sep='\t', index=False)

    # Add metadata
    adata.uns['downsampling'] = {
        'target_depth': args.target_depth,
        'seed': args.seed,
        'source_file': args.input
    }

    # Save
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    print(f"\nSaving to {args.output}...")
    adata.write(args.output)
    print("Done!")


if __name__ == "__main__":
    main()
