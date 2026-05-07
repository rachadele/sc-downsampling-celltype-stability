# Downsampling Experiment

Determine the minimum sequencing depth required for stable cell type calling in Smart-seq v4 single-cell RNA-seq data.

## Goal

At what sequencing depth do cell type assignments stabilize? Is 10x-equivalent depth (~50-100k reads/cell) sufficient for accurate cell typing on Smart-seq data?

## Experimental Design

1. **Download** human Smart-seq v4 reference data from [CellxGene](https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f) via direct curl
2. **Downsample** to various depths using multinomial sampling
3. **Create paired replicates** at each depth (different random seeds)
4. **Run label transfer** using both Seurat and scVI methods
5. **Compare** cell type assignments between replicates
6. **Identify** the minimum depth where agreement stabilizes

### Depth Levels Tested

| Depth | Context |
|-------|---------|
| 1,000 | Very low |
| 5,000 | Low |
| 10,000 | Below typical 10x |
| 15,000 | Low-mid |
| 25,000 | Mid |
| 50,000 | Typical 10x range |
| 100,000 | High 10x |

## Project Structure

```
downsampling-experiment/
├── main.nf              # Nextflow pipeline
├── nextflow.config      # Configuration
├── README.md
├── bin/
│   ├── download_data.py       # Fetch data from CellxGene
│   ├── split_data.py          # Split into reference/query
│   ├── setup.py               # Setup scVI model
│   ├── downsample.py          # Multinomial downsampling
│   ├── predict_scvi.py        # scVI label transfer
│   ├── compare_replicates.py  # Calculate agreement metrics
│   ├── plot_stability.py      # Generate stability curves
│   ├── utils.py               # Python helper functions
│   ├── predict_seurat.R       # (Disabled) Seurat label transfer
│   ├── seurat_functions.R     # (Disabled) R helper functions
│   └── seurat_processing.R    # (Disabled) R processing
├── data/                # Downloaded/input data
└── results/             # Pipeline output
```

## Requirements

### Conda Environments

- **scanpyenv**: Python environment with scanpy, scvi-tools, cellxgene-census
- **r4.3**: R environment with Seurat, SeuratDisk

### Software

- Nextflow >= 21.04.0

## Usage

### Run Full Pipeline

```bash
cd /space/grp/rschwartz/rschwartz/downsampling-experiment
nextflow run main.nf
```

### With Existing Data

```bash
nextflow run main.nf --input_data /path/to/your/smartseq.h5ad
```

### Resume After Interruption

```bash
nextflow run main.nf -resume
```

### Custom Parameters

```bash
nextflow run main.nf \
    --depth_levels '[50000, 100000, 500000]' \
    --label_key 'cell_type' \
    --outdir 'my_results'
```

## Pipeline Steps

```
┌─────────────────┐
│  Download Data  │  (or use --input_data)
└────────┬────────┘
         │
    ┌────┴────┐
    ▼         ▼
┌───────┐ ┌───────┐
│Seurat │ │ scVI  │  Prepare references
│  ref  │ │  ref  │
└───┬───┘ └───┬───┘
    │         │
    └────┬────┘
         │
         ▼
┌─────────────────┐
│   Downsample    │  For each depth × 2 seeds
└────────┬────────┘
         │
    ┌────┴────┐
    ▼         ▼
┌───────┐ ┌───────┐
│Seurat │ │ scVI  │  Label transfer
│predict│ │predict│
└───┬───┘ └───┬───┘
    │         │
    └────┬────┘
         │
         ▼
┌─────────────────┐
│    Compare      │  Agreement between rep1 vs rep2
│   Replicates    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Plot Results   │  Stability curves
└─────────────────┘
```

## Output

### Main Results (`results/results/`)

- `agreement_curves.png` - Agreement vs depth for both methods
- `per_type_heatmap_seurat.png` - Per-cell-type stability (Seurat)
- `per_type_heatmap_scvi.png` - Per-cell-type stability (scVI)
- `stability_summary.tsv` - Minimum depth for 95% agreement
- `all_metrics.tsv` - All comparison metrics

### Intermediate Files

- `results/downsampled/` - Downsampled h5ad files
- `results/predictions/` - Cell type predictions per method/depth/replicate
- `results/comparisons/` - Pairwise comparison metrics

## Methods

### Downsampling

Multinomial sampling treats each cell's gene expression as a probability distribution and samples N reads:

```python
probs = cell_counts / total_counts
sampled = np.random.multinomial(target_depth, probs)
```

### Label Transfer

**Seurat**: `FindTransferAnchors()` + `TransferData()` with SCTransform normalization

**scVI**: Train RandomForest classifier on scVI latent embeddings from reference, predict on query

### Agreement Metrics

- **Percent Agreement**: Fraction of cells with identical assignments
- **Adjusted Rand Index (ARI)**: Agreement corrected for chance
- **Cohen's Kappa**: Inter-rater reliability measure

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `depth_levels` | [10k, 50k, 100k, 250k, 500k, 1M] | Reads per cell to test |
| `seeds` | [42, 123] | Random seeds for paired replicates |
| `label_key` | "cell_type" | Column with cell type annotations |
| `stability_threshold` | 0.95 | Agreement threshold for "stable" |
| `dims` | 50 | PCA dimensions for Seurat |

## Data Source

Direct download URL:
```
https://datasets.cellxgene.cziscience.com/dfba29ea-d368-44c2-bb35-69d3e8f730ce.h5ad
```

Collection: [Human Multiple Cortical Areas SMART-seq](https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f)

## References

- Label transfer methods adapted from [eval-references](../eval-references/)
