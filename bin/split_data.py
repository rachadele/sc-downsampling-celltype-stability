#!/usr/bin/env python3
import argparse
import os
import scanpy as sc
import numpy as np
from sklearn.model_selection import train_test_split

def parse_arguments():
    parser = argparse.ArgumentParser(description="Split AnnData into reference and query sets")
    parser.add_argument('--input', type=str, required=False, default='/space/grp/rschwartz/rschwartz/downsampling-experiment/homo_sapiens/Multiple Cortical Areas SMART-seq/2024-07-01/data/reference.h5ad', help='Input h5ad file')
    parser.add_argument('--ref_output', type=str, help='Output reference h5ad file')
    parser.add_argument('--query_output', type=str, help='Output query h5ad file')
    parser.add_argument('--test_size', type=float, default=0.3, help='Fraction of cells for query set')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--stratify_key', type=str, default="cell_type", help='Key in obs to stratify the split')
    args, unknown = parser.parse_known_args()
    return args
def main():
    args = parse_arguments()
    
    print(f"Reading {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    print(f"Splitting data with test_size={args.test_size}...")
    
    indices = np.arange(adata.n_obs)
               
    train_idx, test_idx = train_test_split(
        indices, 
        test_size=args.test_size, 
        random_state=args.seed,
        stratify=adata.obs[args.stratify_key]  # stratifying is too annoying with rare cell types
    )
    
    adata_ref = adata[train_idx].copy()
    adata_query = adata[test_idx].copy()
    
    print(f"Reference set: {adata_ref.n_obs} cells")
    print(f"Query set: {adata_query.n_obs} cells")
    
    print(f"Saving reference to {args.ref_output}...")
    adata_ref.write_h5ad(args.ref_output)
    
    print(f"Saving query to {args.query_output}...")
    adata_query.write_h5ad(args.query_output)
    
    print("Done!")

if __name__ == "__main__":
    main()
