#!/usr/bin/env python3
"""
Compare cell type assignments between two downsampled replicates.

Computes agreement metrics to assess stability of cell type calling
at a given sequencing depth.
"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.metrics import (
    f1_score,
    adjusted_rand_score,
    confusion_matrix,
    cohen_kappa_score
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compare cell type assignments between two replicates"
    )
    parser.add_argument(
        '--rep1',
        type=str,
        required=True,
        help='Predictions file for replicate 1 (TSV with predicted_label column)'
    )
    parser.add_argument(
        '--rep2',
        type=str,
        required=True,
        help='Predictions file for replicate 2 (TSV with predicted_label column)'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output file for comparison metrics (TSV)'
    )
    parser.add_argument(
        '--depth',
        type=int,
        required=True,
        help='Sequencing depth (counts per cell) for this comparison'
    )
    parser.add_argument(
        '--method',
        type=str,
        required=True,
        help='Method name (seurat or scvi)'
    )
    parser.add_argument(
        '--confusion_matrix',
        type=str,
        default=None,
        help='Optional: output path for confusion matrix'
    )
    parser.add_argument(
        '--seed1',
        type=int,
        required=True,
        help='Random seed used for replicate 1'
    )
    parser.add_argument(
        '--seed2',
        type=int,
        required=True,
        help='Random seed used for replicate 2'
    )
    return parser.parse_args()


def compute_agreement_metrics(labels1, labels2, all_labels=None):
    """
    Compute various agreement metrics between two sets of labels.

    Parameters
    ----------
    labels1 : array-like
        Cell type labels from replicate 1
    labels2 : array-like
        Cell type labels from replicate 2
    all_labels : array-like, optional
        All possible labels to compute per-type metrics for (ensures consistency)

    Returns
    -------
    metrics : dict
        Dictionary of agreement metrics
    """
    labels1 = np.array(labels1)
    labels2 = np.array(labels2)

    # Filter out any cells where either label is missing
    valid_mask = ~(pd.isna(labels1) | pd.isna(labels2))
    labels1 = labels1[valid_mask]
    labels2 = labels2[valid_mask]

    n_cells = len(labels1)

    # Percent agreement (exact match)
    exact_matches = (labels1 == labels2).sum()
    percent_agreement = exact_matches / n_cells if n_cells > 0 else 0

    # Adjusted Rand Index
    ari = adjusted_rand_score(labels1, labels2)

    # Cohen's Kappa
    kappa = cohen_kappa_score(labels1, labels2)

    # Per-cell-type accuracy and counts
    # Use provided all_labels if given, otherwise compute from data
    if all_labels is None:
        unique_labels = np.unique(np.concatenate([labels1, labels2]))
    else:
        unique_labels = all_labels

    per_type_f1 = {}
    per_type_counts = {}

    for label in unique_labels:
        # Binary labels: 1 if cell is this type, 0 otherwise
        binary1 = (labels1 == label).astype(int)
        binary2 = (labels2 == label).astype(int)

        # Cells called as this type in rep1 and rep2
        n_rep1 = binary1.sum()
        n_rep2 = binary2.sum()
        per_type_counts[label] = {'rep1': int(n_rep1), 'rep2': int(n_rep2)}

        # Compute F1 score - will be 0 if either has no predictions for this type
        # zero_division=0 ensures F1=0 when there are no positive predictions
        per_type_f1[label] = f1_score(binary1, binary2, zero_division=0)

    return {
        'n_cells': n_cells,
        'percent_agreement': percent_agreement,
        'adjusted_rand_index': ari,
        'cohen_kappa': kappa,
        'n_cell_types': len(unique_labels),
        'per_type_f1': per_type_f1,
        'per_type_counts': per_type_counts
    }


def main():
    args = parse_arguments()

    # Load predictions

    rep1_df = pd.read_csv(args.rep1, sep='\t')
    rep2_df = pd.read_csv(args.rep2, sep='\t')

    # Verify same cells (by index order - assumes same cell order)
    if len(rep1_df) != len(rep2_df):
        raise ValueError(f"Replicate files have different number of cells: {len(rep1_df)} vs {len(rep2_df)}")

    # Get predictions and original labels
    labels1 = rep1_df['predicted_label'].values
    labels2 = rep2_df['predicted_label'].values
    orig1 = rep1_df['original_label'].values if 'original_label' in rep1_df.columns else None
    orig2 = rep2_df['original_label'].values if 'original_label' in rep2_df.columns else None

    # Collect all unique labels across all comparisons for consistency
    all_labels_list = [labels1, labels2]
    if orig1 is not None:
        all_labels_list.append(orig1)
    if orig2 is not None:
        all_labels_list.append(orig2)
    all_unique_labels = np.unique(np.concatenate(all_labels_list))

    # Compute replicate-to-replicate metrics
    metrics = compute_agreement_metrics(labels1, labels2, all_labels=all_unique_labels)

    # Compute agreement to original labels if available
    orig1_metrics = compute_agreement_metrics(labels1, orig1, all_labels=all_unique_labels) if orig1 is not None else None
    orig2_metrics = compute_agreement_metrics(labels2, orig2, all_labels=all_unique_labels) if orig2 is not None else None



    # Create output dataframe
    output_df = pd.DataFrame([{
        'depth': args.depth,
        'method': args.method,
        'seed1': args.seed1,
        'seed2': args.seed2,
        'n_cells': metrics['n_cells'],
        'percent_agreement': metrics['percent_agreement'],
        'adjusted_rand_index': metrics['adjusted_rand_index'],
        'cohen_kappa': metrics['cohen_kappa'],
        'n_cell_types': metrics['n_cell_types'],
        'rep1_vs_original_percent_agreement': orig1_metrics['percent_agreement'] if orig1_metrics is not None else None,
        'rep1_vs_original_ari': orig1_metrics['adjusted_rand_index'] if orig1_metrics is not None else None,
        'rep2_vs_original_percent_agreement': orig2_metrics['percent_agreement'] if orig2_metrics is not None else None,
        'rep2_vs_original_ari': orig2_metrics['adjusted_rand_index'] if orig2_metrics is not None else None
    }])

    # Save main metrics
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    output_df.to_csv(args.output, sep='\t', index=False)


    # Save per-cell-type F1, including original label comparisons if available
    per_type_output = args.output.replace('.tsv', '_per_type.tsv')
    per_type_rows = [
        {
            'cell_type': ct,
            'f1': f1,
            'count_rep1': metrics['per_type_counts'][ct]['rep1'],
            'count_rep2': metrics['per_type_counts'][ct]['rep2'],
            'depth': args.depth,
            'method': args.method,
            'seed1': args.seed1,
            'seed2': args.seed2,
            'comparison': 'rep1_vs_rep2'
        }
        for ct, f1 in metrics['per_type_f1'].items()
    ]
    if orig1_metrics is not None:
        per_type_rows.extend([
            {
                'cell_type': ct,
                'f1': f1,
                'count_rep1': orig1_metrics['per_type_counts'][ct]['rep1'],
                'count_rep2': orig1_metrics['per_type_counts'][ct]['rep2'],
                'depth': args.depth,
                'method': args.method,
                'seed1': args.seed1,
                'seed2': args.seed2,
                'comparison': 'rep1_vs_original'
            }
            for ct, f1 in orig1_metrics['per_type_f1'].items()
        ])
    if orig2_metrics is not None:
        per_type_rows.extend([
            {
                'cell_type': ct,
                'f1': f1,
                'count_rep1': orig2_metrics['per_type_counts'][ct]['rep1'],
                'count_rep2': orig2_metrics['per_type_counts'][ct]['rep2'],
                'depth': args.depth,
                'method': args.method,
                'seed1': args.seed1,
                'seed2': args.seed2,
                'comparison': 'rep2_vs_original'
            }
            for ct, f1 in orig2_metrics['per_type_f1'].items()
        ])
    per_type_df = pd.DataFrame(per_type_rows)
    per_type_df.to_csv(per_type_output, sep='\t', index=False)


    # Optionally save confusion matrix
    if args.confusion_matrix:
        all_labels = np.unique(np.concatenate([labels1, labels2]))
        cm = confusion_matrix(labels1, labels2, labels=all_labels)
        cm_df = pd.DataFrame(cm, index=all_labels, columns=all_labels)
        cm_df.to_csv(args.confusion_matrix, sep='\t')

        # Also save confusion matrices for original vs replicate predictions if available
        if 'original_label' in rep1_df.columns:
            orig1 = rep1_df['original_label'].values
            all_labels1 = np.unique(np.concatenate([orig1, labels1]))
            cm1 = confusion_matrix(orig1, labels1, labels=all_labels1)
            cm1_df = pd.DataFrame(cm1, index=all_labels1, columns=all_labels1)
            out1 = args.confusion_matrix.replace('.tsv', '_rep1_vs_original.tsv')
            cm1_df.to_csv(out1, sep='\t')
        if 'original_label' in rep2_df.columns:
            orig2 = rep2_df['original_label'].values
            all_labels2 = np.unique(np.concatenate([orig2, labels2]))
            cm2 = confusion_matrix(orig2, labels2, labels=all_labels2)
            cm2_df = pd.DataFrame(cm2, index=all_labels2, columns=all_labels2)
            out2 = args.confusion_matrix.replace('.tsv', '_rep2_vs_original.tsv')
            cm2_df.to_csv(out2, sep='\t')




if __name__ == "__main__":
    main()
