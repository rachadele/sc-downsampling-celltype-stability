#!/usr/bin/env python3
"""
Extract compute times from Nextflow trace file and plot by method and depth.
"""

import argparse
import os
import re
import pandas as pd
import matplotlib.pyplot as plt


def parse_duration(duration_str):
    """Parse Nextflow duration string to seconds."""
    if pd.isna(duration_str):
        return None

    duration_str = str(duration_str).strip()
    total_seconds = 0

    # Match patterns like "1m 30s", "45.5s", "2h 30m", etc.
    hours = re.search(r'(\d+)h', duration_str)
    minutes = re.search(r'(\d+)m', duration_str)
    seconds = re.search(r'([\d.]+)s', duration_str)

    if hours:
        total_seconds += int(hours.group(1)) * 3600
    if minutes:
        total_seconds += int(minutes.group(1)) * 60
    if seconds:
        total_seconds += float(seconds.group(1))

    return total_seconds if total_seconds > 0 else None


def extract_depth_from_workdir(work_dir, hash_str):
    """Extract depth from the .command.sh file in the work directory."""
    # Construct the work directory path
    short_hash = hash_str.split('/')[0][:2]
    long_hash = hash_str.split('/')[1] if '/' in hash_str else hash_str

    cmd_file = os.path.join(work_dir, short_hash, long_hash + '*', '.command.sh')

    # Use glob to find the actual directory
    import glob
    matches = glob.glob(os.path.join(work_dir, short_hash, long_hash + '*', '.command.sh'))

    if not matches:
        return None

    try:
        with open(matches[0], 'r') as f:
            content = f.read()

        # Look for depth in different patterns
        # For PREDICT_SCVI/SEURAT, look for the input file pattern
        depth_match = re.search(r'downsampled_(\d+)_\d+', content)
        if depth_match:
            return int(depth_match.group(1))

        # For DOWNSAMPLE, look for --target_depth
        depth_match = re.search(r'--target_depth\s+(\d+)', content)
        if depth_match:
            return int(depth_match.group(1))

    except Exception as e:
        pass

    return None


def main():
    parser = argparse.ArgumentParser(description='Plot compute times from Nextflow trace')
    parser.add_argument('--trace', type=str, default="/space/grp/rschwartz/rschwartz/downsampling-experiment/homo_sapiens/Multiple Cortical Areas SMART-seq/2024-07-01/trace.txt", help='Path to trace.txt file')
    parser.add_argument('--work_dir', type=str, default="/space/grp/rschwartz/rschwartz/downsampling-experiment/work", help='Path to Nextflow work directory')
    parser.add_argument('--output', type=str, default="/space/grp/rschwartz/rschwartz/downsampling-experiment/homo_sapiens/Multiple Cortical Areas SMART-seq/2024-07-01/results/compute_times.png", help='Output plot file')
    parser.add_argument('--output_tsv', type=str, default="/space/grp/rschwartz/rschwartz/downsampling-experiment/homo_sapiens/Multiple Cortical Areas SMART-seq/2024-07-01/results/compute_times.tsv", help='Output TSV file with timing data')
    args = parser.parse_args()

    # Read trace file
    df = pd.read_csv(args.trace, sep='\t')

    # Filter to PREDICT_SCVI and PREDICT_SEURAT tasks
    predict_df = df[df['name'].str.contains('PREDICT_SCVI|PREDICT_SEURAT', regex=True)].copy()

    # Parse realtime to seconds
    predict_df['realtime_seconds'] = predict_df['realtime'].apply(parse_duration)

    # Extract method from name
    predict_df['method'] = predict_df['name'].apply(
        lambda x: 'scvi' if 'SCVI' in x else 'seurat'
    )

    # Extract depth from work directories
    predict_df['depth'] = predict_df['hash'].apply(
        lambda h: extract_depth_from_workdir(args.work_dir, h)
    )

    # Remove rows where we couldn't get depth
    predict_df = predict_df.dropna(subset=['depth', 'realtime_seconds'])

    # Save TSV if requested
    if args.output_tsv:
        output_cols = ['name', 'method', 'depth', 'realtime_seconds', 'realtime', 'peak_rss']
        predict_df[output_cols].to_csv(args.output_tsv, sep='\t', index=False)

    # Aggregate by method and depth (average over seeds)
    agg_df = predict_df.groupby(['method', 'depth']).agg({
        'realtime_seconds': 'mean'
    }).reset_index()

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {'scvi': 'tab:blue', 'seurat': 'tab:red'}
    markers = {'scvi': 'o', 'seurat': 's'}

    for method in ['scvi', 'seurat']:
        method_df = agg_df[agg_df['method'] == method].sort_values('depth')
        ax.plot(
            method_df['depth'],
            method_df['realtime_seconds'],
            marker=markers[method],
            color=colors[method],
            label=method.upper(),
            linewidth=2,
            markersize=8
        )

    ax.set_xlabel('Reads per Cell', fontsize=12)
    ax.set_ylabel('Compute Time (seconds)', fontsize=12)
    ax.set_title('Prediction Compute Time vs Sequencing Depth', fontsize=14)
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {args.output}")


if __name__ == '__main__':
    main()
