#!/usr/bin/env python3
"""
scVI-based label transfer for downsampling experiment.

Uses scVI embeddings + RandomForest classifier to transfer cell type
labels from reference to downsampled query.
"""

import argparse
import os
import random
import warnings
warnings.filterwarnings("ignore")
from utils import process_query
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
from sklearn.ensemble import RandomForestClassifier


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="scVI label transfer for downsampling experiment"
    )
    parser.add_argument(
            '--batch_components',
            type=str,
            nargs='+',
            default=['assay', 'tissue', 'dataset_id'],
            help='List of obs columns to use for batch construction (space-separated)'
        )
    parser.add_argument(
        '--ref_path',
        type=str,
        default='/space/grp/rschwartz/rschwartz/downsampling-experiment/work/1d/1cb475dd33e923befe7dce53edc722/reference.h5ad',
        help='Path to reference h5ad file (with scVI embeddings)'
    )
    parser.add_argument(
        '--query_path',
        type=str,
        default='/space/grp/rschwartz/rschwartz/downsampling-experiment/work/1d/1cb475dd33e923befe7dce53edc722/downsampled_100000_42.h5ad',
        help='Path to downsampled query h5ad file'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='/space/grp/rschwartz/rschwartz/downsampling-experiment/work/1d/1cb475dd33e923befe7dce53edc722',
        help='Output directory'
    )
    parser.add_argument(
        '--label_key',
        type=str,
        default='cell_type',
        help='Column name for cell type labels'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed'
    )
    parser.add_argument(
        '--model_path',
        type=str,
        default='/space/grp/rschwartz/rschwartz/downsampling-experiment/work/1d/1cb475dd33e923befe7dce53edc722/scvi-homo_sapiens-2024-07-01',
        help='Path to pre-trained scVI model directory or file'
    )
    parser.add_argument(
        '--batch_key',
        type=str,
        default='batch',
        help='Batch key for scVI model'
    )
    args, _ = parser.parse_known_args()
    return args



def train_and_predict(ref_embeddings, ref_labels, query_embeddings, seed=42):
    """
    Train RandomForest on reference embeddings and predict query labels.

    Parameters
    ----------
    ref_embeddings : np.array
        Reference scVI embeddings (n_ref_cells x n_latent)
    ref_labels : np.array
        Reference cell type labels
    query_embeddings : np.array
        Query scVI embeddings (n_query_cells x n_latent)
    seed : int
        Random seed

    Returns
    -------
    predictions : np.array
        Predicted labels for query cells
    probabilities : np.array
        Prediction probabilities (n_query_cells x n_classes)
    class_labels : np.array
        Class label names
    """
    print("Training RandomForest classifier...")

    rfc = RandomForestClassifier(
        n_estimators=100,
        class_weight='balanced',
        random_state=seed,
        n_jobs=-1
    )

    rfc.fit(ref_embeddings, ref_labels)

    print("Predicting query labels...")
    predictions = rfc.predict(query_embeddings)
    probabilities = rfc.predict_proba(query_embeddings)
    class_labels = rfc.classes_

    return predictions, probabilities, class_labels


def main():
    args = parse_arguments()

    # Set seeds
    random.seed(args.seed)
    np.random.seed(args.seed)
    scvi.settings.seed = args.seed

    os.makedirs(args.output_dir, exist_ok=True)

    # Load reference
    print(f"Loading reference: {args.ref_path}")
    ref = ad.read_h5ad(args.ref_path, backed='r')

    # Check for scVI embeddings in reference
    if 'scvi' not in ref.obsm:
        raise ValueError("Reference must have scVI embeddings in .obsm['scvi']")

    ref_embeddings = ref.obsm['scvi']
    ref_labels = ref.obs[args.label_key].values

    print(f"Reference: {ref.n_obs} cells, {len(np.unique(ref_labels))} cell types")

    # Load query
    print(f"Loading query: {args.query_path}")
    query = ad.read_h5ad(args.query_path)

    print(f"Query: {query.n_obs} cells")

    # Always compute scVI embeddings for the query using the reference model
    print("Computing new scVI embeddings for query using reference model...")
    
    # construct batch from user-specified obs columns
    batch_components = args.batch_components
    query.obs['batch'] = query.obs[batch_components].astype(str).agg('_'.join, axis=1)
    
    query = process_query(query, args.model_path, batch_key=args.batch_key, seed=args.seed)
     
    query_embeddings = query.obsm['scvi']
    # Train and predict
    predictions, probabilities, class_labels = train_and_predict(
        ref_embeddings, ref_labels, query_embeddings, args.seed
    )

    # Create output dataframe
    output_df = pd.DataFrame({
        'cell_id': query.obs.index,
        'predicted_label': predictions,
        'max_score': probabilities.max(axis=1)
    })

    # Add original labels
    output_df['original_label'] = query.obs[args.label_key].values

    # Create probability dataframe
    prob_df = pd.DataFrame(probabilities, columns=class_labels)

    # Save outputs
    query_name = os.path.basename(args.query_path).replace('.h5ad', '')

    pred_path = os.path.join(args.output_dir, f"{query_name}_scvi_predictions.tsv")
    output_df.to_csv(pred_path, sep='\t', index=False)
    print(f"Saved predictions to: {pred_path}")

    scores_path = os.path.join(args.output_dir, f"{query_name}_scvi_scores.tsv")
    prob_df.to_csv(scores_path, sep='\t', index=False)
    print(f"Saved scores to: {scores_path}")

    print("Done!")


if __name__ == "__main__":
    main()
