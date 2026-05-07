#!/usr/bin/env python3
"""
Download human Smart-seq data from CellxGene Census for downsampling experiment.

"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import cellxgene_census
import scanpy as sc
import numpy as np
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Download data from CellxGene Census"
    )
    parser.add_argument(
        '--dataset_name',
        type=str,
        default="Multiple Cortical Areas SMART-seq",
        help='Dataset name to search for in CellxGene Census'
    )
    parser.add_argument(
        '--organism',
        type=str,
        default='homo_sapiens',
        help='Organism (homo_sapiens or mus_musculus)'
    )
    parser.add_argument(
        '--census_version',
        type=str,
        default='2024-07-01',
        help='Census version'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='output.h5ad',
        help='Output file path'
    )
    parser.add_argument(
        '--n_cells',
        type=int,
        default=5000,
        help='Number of cells to subsample (default: all cells)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed'
    )
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        


def main():
    args = parse_arguments()
    np.random.seed(args.seed)

    # Map organism argument to census format
    organism_map = {
        'homo_sapiens': 'Homo sapiens',
        'mus_musculus': 'Mus musculus'
    }
    organism_display = organism_map.get(args.organism, args.organism)

    print(f"Opening CellxGene Census version {args.census_version}...")
    census = cellxgene_census.open_soma(census_version=args.census_version)

    # Get dataset info
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()

    # Filter by dataset name
    matching_datasets = dataset_info[
        dataset_info['dataset_title'].str.contains(args.dataset_name, case=False, na=False)
    ]

    if matching_datasets.empty:
        raise ValueError(f"No datasets found matching '{args.dataset_name}'") 

    print(f"Found {len(matching_datasets)} matching dataset(s):")
    print(matching_datasets[['dataset_id', 'dataset_title', 'collection_name']].to_string())

    # Get the dataset IDs
    dataset_ids = matching_datasets['dataset_id'].tolist()

    # Query cells from these datasets
    print("\nQuerying cells from census...")

    obs_df = census["census_data"][args.organism].obs.read(
        value_filter=f"dataset_id in {dataset_ids} and is_primary_data == True",
        column_names=[
            "soma_joinid", "cell_type", "assay", "tissue",
            "dataset_id", "donor_id", "suspension_type",
            "disease", "development_stage"
        ]
    ).concat().to_pandas()

    print(f"Found {len(obs_df)} cells")

    # Subsample cells if requested
    if args.n_cells:
        print(f"\nSubsampling to {args.n_cells} cells...")
        obs_df = obs_df.sample(n=args.n_cells, random_state=args.seed)
        print(f"Subsampled to {len(obs_df)} cells")

    print(f"\nCell type distribution:")
    print(obs_df['cell_type'].value_counts())

    # Get the expression data with scVI embeddings
    print("\nFetching expression data (this may take a while)...")

    cell_ids = obs_df['soma_joinid'].tolist()

    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism_display,
        obs_coords=cell_ids,
        obs_embeddings=["scvi"],
        column_names={
            "obs": ["soma_joinid", "cell_type", "assay", "tissue",
                   "dataset_id", "donor_id", "disease", "development_stage"],
            "var": ["feature_id", "feature_name"]
        }
    )

    # merge dataset info on soma_joinid so that dataset_title is included
    
    adata.obs = adata.obs.merge(
        matching_datasets,
        on='dataset_id',
        how='left'
    )


    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Calculate total counts per cell for reference
    if hasattr(adata.X, 'toarray'):
        adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    else:
        adata.obs['total_counts'] = adata.X.sum(axis=1)

    print(f"\nFinal dataset shape: {adata.shape}")
    print(f"Total counts range: {adata.obs['total_counts'].min():.0f} - {adata.obs['total_counts'].max():.0f}")
    print(f"Mean counts per cell: {adata.obs['total_counts'].mean():.0f}")
    print(f"Median counts per cell: {adata.obs['total_counts'].median():.0f}")

    # Save the data
    print(f"\nSaving to {args.output}...")
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    adata.write(args.output)

    print("Done!")
    census.close()

    return adata


if __name__ == "__main__":
    main()
