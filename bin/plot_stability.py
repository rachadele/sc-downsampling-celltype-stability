#!/usr/bin/env python3
"""
Plot stability curves for the downsampling experiment.

Aggregates comparison metrics across all depth levels and creates
visualizations to identify the minimum depth for stable cell typing.
"""

import argparse
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot stability curves from downsampling experiment"
    )
    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Directory containing comparison metric files'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory for plots'
    )
    parser.add_argument(
        '--stability_threshold',
        type=float,
        default=0.95,
        help='Threshold for "stable" agreement (default: 0.95)'
    )
    return parser.parse_args()


def find_stability_point(df, metric='percent_agreement', threshold=0.95):
    """
    Find the minimum depth where metric exceeds threshold.

    Returns the depth and metric value, or None if threshold never reached.
    """
    df_sorted = df.sort_values('depth')

    for _, row in df_sorted.iterrows():
        if row[metric] >= threshold:
            return row['depth'], row[metric]

    return None, None


def plot_agreement_curves(df, output_path, threshold=0.95):
    """
    Plot ARI vs depth for both methods.
    Shows stability (rep1 vs rep2) and accuracy (vs original).
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    # Define colors for each method
    method_colors = {'scvi': 'tab:blue', 'seurat': 'tab:red'}
    method_markers = {'scvi': 'o', 'seurat': 's'}

    plot_configs = [
        ('adjusted_rand_index', 'ARI: Rep1 vs Rep2 (Stability)'),
        ('rep1_vs_original_ari', 'ARI: Rep1 vs Original (Accuracy)'),
        ('rep2_vs_original_ari', 'ARI: Rep2 vs Original (Accuracy)')
    ]

    for idx, (ax, (metric, title)) in enumerate(zip(axes, plot_configs)):
        if metric not in df.columns:
            ax.set_visible(False)
            continue

        for method in df['method'].unique():
            method_df = df[df['method'] == method].sort_values('depth')
            color = method_colors.get(method, 'tab:gray')
            marker = method_markers.get(method, 'o')

            ax.plot(
                method_df['depth'],
                method_df[metric],
                marker=marker,
                color=color,
                label=f'{method.upper()}',
                linewidth=2,
                markersize=8
            )

        ax.set_xlabel('Reads per Cell', fontsize=12)
        if idx == 0:
            ax.set_ylabel('Adjusted Rand Index', fontsize=12)
        ax.set_title(title, fontsize=14)
        ax.set_xscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)

    axes[0].set_ylim(-1, 1.05)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved agreement curves to: {output_path}")


def plot_per_type_heatmap(df, output_path):
    """
    Create stacked heatmaps showing per-cell-type agreement across depths.
    Shows rep1 vs rep2, rep1 vs original, and rep2 vs original.
    """
    # Determine which column to use (f1 or accuracy for backwards compatibility)
    value_col = 'f1' if 'f1' in df.columns else 'accuracy'
    metric_label = 'F1' if value_col == 'f1' else 'Accuracy'

    comparisons = ['rep1_vs_rep2', 'rep1_vs_original', 'rep2_vs_original']
    comparison_titles = {
        'rep1_vs_rep2': f'{metric_label}: Rep1 vs Rep2',
        'rep1_vs_original': f'{metric_label}: Rep1 vs Original',
        'rep2_vs_original': f'{metric_label}: Rep2 vs Original'
    }

    for method in df['method'].unique():
        method_df = df[df['method'] == method]

        # Check which comparisons are available
        available_comparisons = [c for c in comparisons if c in method_df['comparison'].unique()]
        n_comparisons = len(available_comparisons)

        if n_comparisons == 0:
            continue

        # Get consistent cell type ordering based on rep1_vs_rep2 mean score
        rep1_rep2_df = method_df[method_df['comparison'] == 'rep1_vs_rep2']
        pivot_for_order = rep1_rep2_df.pivot(index='cell_type', columns='depth', values=value_col)
        cell_type_order = pivot_for_order.mean(axis=1).sort_values(ascending=False).index.tolist()

        # Stack horizontally (1 row, n_comparisons columns)
        n_cell_types = len(cell_type_order)
        fig, axes = plt.subplots(1, n_comparisons, figsize=(6 * n_comparisons, max(8, n_cell_types * 0.4)))
        if n_comparisons == 1:
            axes = [axes]

        for idx, (ax, comparison) in enumerate(zip(axes, available_comparisons)):
            comp_df = method_df[method_df['comparison'] == comparison]

            # Pivot to create heatmap data
            pivot_df = comp_df.pivot(
                index='cell_type',
                columns='depth',
                values=value_col
            )

            # Reorder to match consistent ordering, fill missing with 0
            pivot_df = pivot_df.reindex(cell_type_order).fillna(0)

            # Only show colorbar on the last plot
            cbar = (idx == n_comparisons - 1)

            # Only show y-tick labels on the first plot
            yticklabels = True if idx == 0 else False

            sns.heatmap(
                pivot_df,
                annot=True,
                fmt='.2f',
                cmap='RdYlGn',
                vmin=0,
                vmax=1,
                ax=ax,
                cbar=cbar,
                cbar_kws={'label': metric_label} if cbar else {},
                yticklabels=yticklabels
            )

            # Only show y-axis label on the first plot
            if idx == 0:
                ax.set_ylabel('Cell Type', fontsize=12)
            else:
                ax.set_ylabel('')

            ax.set_xlabel('Reads per Cell', fontsize=12)
            ax.set_title(f'{comparison_titles[comparison]}', fontsize=12)

        plt.suptitle(f'Per-Cell-Type {metric_label} ({method.upper()})', fontsize=14)
        plt.tight_layout()
        method_output = output_path.replace('.png', f'_{method}.png')
        plt.savefig(method_output, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved per-type heatmap to: {method_output}")


def save_per_type_counts(df, output_path):
    """
    Save a table showing cell type counts for each method and depth.
    """
    # Filter to rep1_vs_rep2 comparisons
    df_filtered = df[df['comparison'] == 'rep1_vs_rep2'].copy()

    # Select relevant columns (use f1 or accuracy depending on what's available)
    value_col = 'f1' if 'f1' in df_filtered.columns else 'accuracy'
    cols = ['method', 'depth', 'cell_type', value_col]
    if 'count_rep1' in df_filtered.columns:
        cols.extend(['count_rep1', 'count_rep2'])

    counts_df = df_filtered[cols].copy()
    counts_df = counts_df.sort_values(['method', 'depth', 'cell_type'])
    counts_df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved per-type counts to: {output_path}")

    # Print summary of cell type counts by method
    if 'count_rep1' in counts_df.columns:
        print(f"\nCell type counts by method:")
        for method in counts_df['method'].unique():
            method_df = counts_df[counts_df['method'] == method]
            print(f"  {method.upper()}:")
            # Show counts for the first depth as example
            first_depth = method_df['depth'].min()
            depth_df = method_df[method_df['depth'] == first_depth]
            for _, row in depth_df.iterrows():
                print(f"    {row['cell_type']}: rep1={row['count_rep1']}, rep2={row['count_rep2']}")

    return counts_df


def create_summary_table(df, threshold=0.95):
    """
    Create summary table with stability points for each method using ARI.
    """
    summary = []

    for method in df['method'].unique():
        method_df = df[df['method'] == method]

        # Find stability point using ARI
        depth, _ = find_stability_point(
            method_df, 'adjusted_rand_index', threshold
        )

        max_ari = method_df['adjusted_rand_index'].max()
        min_depth_tested = method_df['depth'].min()
        max_depth_tested = method_df['depth'].max()

        # Get accuracy metrics if available
        max_rep1_orig_ari = method_df['rep1_vs_original_ari'].max() if 'rep1_vs_original_ari' in method_df.columns else None
        max_rep2_orig_ari = method_df['rep2_vs_original_ari'].max() if 'rep2_vs_original_ari' in method_df.columns else None

        # Get seeds used if available
        seed1 = method_df['seed1'].iloc[0] if 'seed1' in method_df.columns else None
        seed2 = method_df['seed2'].iloc[0] if 'seed2' in method_df.columns else None

        summary.append({
            'method': method,
            'seed1': seed1,
            'seed2': seed2,
            f'min_depth_for_{int(threshold*100)}%_ari': depth,
            'max_stability_ari': max_ari,
            'max_rep1_vs_original_ari': max_rep1_orig_ari,
            'max_rep2_vs_original_ari': max_rep2_orig_ari,
            'min_depth_tested': min_depth_tested,
            'max_depth_tested': max_depth_tested
        })

    return pd.DataFrame(summary)


def main():
    args = parse_arguments()

    os.makedirs(args.output_dir, exist_ok=True)

    # Find all comparison metric files
    metric_files = []
    per_type_files = []

    for root, dirs, files in os.walk(args.input_dir):
        for f in files:
            if f.endswith('.tsv'):
                full_path = os.path.join(root, f)
                if '_per_type' in f:
                    per_type_files.append(full_path)
                elif 'comparison' in f or 'metrics' in f:
                    metric_files.append(full_path)

    if not metric_files:
        # Try to find files matching pattern
        import glob
        metric_files = glob.glob(os.path.join(args.input_dir, '**/comparison*.tsv'), recursive=True)
        metric_files += glob.glob(os.path.join(args.input_dir, '**/metrics*.tsv'), recursive=True)
        per_type_files = glob.glob(os.path.join(args.input_dir, '**/*_per_type*.tsv'), recursive=True)

    print(f"Found {len(metric_files)} metric files")
    print(f"Found {len(per_type_files)} per-type files")

    if not metric_files:
        print("No metric files found. Please check input directory.")
        return

    # Load and combine all metric files
    dfs = []
    for f in metric_files:
        try:
            df = pd.read_csv(f, sep='\t')
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not load {f}: {e}")

    if not dfs:
        print("No data loaded. Exiting.")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"\nLoaded {len(combined_df)} comparison results")
    print(f"Depths: {sorted(combined_df['depth'].unique())}")
    print(f"Methods: {combined_df['method'].unique()}")

    # Plot agreement curves
    plot_agreement_curves(
        combined_df,
        os.path.join(args.output_dir, 'agreement_curves.png'),
        threshold=args.stability_threshold
    )

    # Load per-type data if available
    if per_type_files:
        per_type_dfs = []
        for f in per_type_files:
            try:
                df = pd.read_csv(f, sep='\t')
                per_type_dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not load {f}: {e}")

        if per_type_dfs:
            per_type_combined = pd.concat(per_type_dfs, ignore_index=True)
            plot_per_type_heatmap(
                per_type_combined,
                os.path.join(args.output_dir, 'per_type_heatmap.png')
            )
            save_per_type_counts(
                per_type_combined,
                os.path.join(args.output_dir, 'per_type_counts.tsv')
            )

    # Create summary table
    summary_df = create_summary_table(combined_df, args.stability_threshold)
    summary_path = os.path.join(args.output_dir, 'stability_summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"\nSaved summary to: {summary_path}")

    print("\n=== Summary ===")
    print(summary_df.to_string(index=False))

    # Save combined data
    combined_path = os.path.join(args.output_dir, 'all_metrics.tsv')
    combined_df.to_csv(combined_path, sep='\t', index=False)
    print(f"\nSaved combined metrics to: {combined_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
