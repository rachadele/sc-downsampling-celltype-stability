# Downsampling Experiment Implementation Plan

## Goal
Determine the minimum sequencing depth required for stable cell type calling in Smart-seq single-cell RNA-seq data.

## Experimental Design

### Key Question
At what read depth do cell type assignments stabilize between two independent downsampling replicates?

### Approach
1. For each depth level X, downsample the original counts **twice** (independently) to create two pseudo-replicates
2. Run both Seurat and scVI-based label transfer on each replicate
3. Compare cell type assignments between the two replicates
4. Find the minimum depth where agreement is high (e.g., >95%)

---

## Implementation Steps

### 1. Data Acquisition

**Script**: `bin/download_data.py`

```
- Use cellxgene_census to download human Smart-seq data
- Collection: https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f
- Ensure raw counts are available (not normalized)
- Save as .h5ad file
```

**Output**: `data/human_smartseq_cells.h5ad`

---

### 2. Downsampling Function

**Script**: `bin/downsample.py`

The downsampling must happen at the **read level**, treating counts as draws from a multinomial distribution:

```python
def downsample_counts(adata, target_reads_per_cell, seed):
    """
    Downsample counts to target depth using multinomial sampling.

    For each cell:
    1. Calculate total UMI/reads in original
    2. Treat gene proportions as multinomial probabilities
    3. Sample `target_reads_per_cell` reads from this distribution
    """
    import numpy as np
    from scipy.sparse import lil_matrix

    np.random.seed(seed)
    downsampled = lil_matrix(adata.X.shape)

    for i in range(adata.n_obs):
        cell_counts = adata.X[i, :].toarray().flatten()
        total_counts = cell_counts.sum()

        if total_counts == 0:
            continue

        # Normalize to get probabilities
        probs = cell_counts / total_counts

        # Sample from multinomial
        n_samples = min(int(target_reads_per_cell), int(total_counts))
        sampled_counts = np.random.multinomial(n_samples, probs)
        downsampled[i, :] = sampled_counts

    return downsampled.tocsr()
```

**Input**: Full-depth h5ad, target depth, random seed
**Output**: Downsampled h5ad

---

### 3. Label Transfer Methods (from eval-references)

**Seurat-based** (`bin/predict_seurat.R`):
- Uses `FindTransferAnchors()` + `TransferData()` from Seurat
- Requires SCTransform normalization
- Already implemented in eval-references

**scVI-based** (`bin/predict_scvi.py`):
- Uses scVI embeddings + RandomForest classifier
- Trains on reference scVI embeddings
- Predicts on query scVI embeddings
- Already implemented in eval-references

For this experiment, we'll use the **non-downsampled** reference and **downsampled** queries.

---

### 4. Nextflow Pipeline

**Script**: `main.nf`

```nextflow
// Depth levels to test (reads per cell)
params.depth_levels = [100, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000]
params.n_replicates = 2  // Create 2 downsampled versions per depth
params.n_cells = 1000
params.ref_key = "cell_type"

workflow {
    // 1. Download and prepare data
    DOWNLOAD_DATA()

    // 2. Subsample to N cells
    SUBSAMPLE_CELLS(DOWNLOAD_DATA.out, params.n_cells)

    // 3. For each depth level, create 2 downsampled replicates
    depth_ch = Channel.from(params.depth_levels)
    replicate_ch = Channel.from(1..params.n_replicates)

    DOWNSAMPLE(SUBSAMPLE_CELLS.out, depth_ch.combine(replicate_ch))
    // Output: tuples of (depth, replicate_id, downsampled_h5ad)

    // 4. Run Seurat label transfer
    SEURAT_PREDICT(DOWNSAMPLE.out)

    // 5. Run scVI label transfer
    SCVI_PREDICT(DOWNSAMPLE.out)

    // 6. Compare replicates at each depth
    COMPARE_REPLICATES(
        SEURAT_PREDICT.out.groupTuple(by: 0),  // group by depth
        SCVI_PREDICT.out.groupTuple(by: 0)
    )

    // 7. Aggregate and plot results
    PLOT_RESULTS(COMPARE_REPLICATES.out.collect())
}
```

---

### 5. Comparison Metrics

**Script**: `bin/compare_replicates.py`

```python
from sklearn.metrics import adjusted_rand_score, accuracy_score
import pandas as pd

def compare_assignments(pred1, pred2):
    """
    Compare cell type assignments between two replicates.

    Returns:
    - percent_agreement: fraction of cells with identical calls
    - adjusted_rand_index: ARI accounting for chance
    - confusion_matrix: which types get confused
    """
    agreement = (pred1 == pred2).mean()
    ari = adjusted_rand_score(pred1, pred2)

    return {
        'percent_agreement': agreement,
        'adjusted_rand_index': ari
    }
```

---

### 6. Analysis and Visualization

**Script**: `bin/plot_stability.py`

Create plots showing:
1. **Agreement vs Depth** - Line plot for Seurat and scVI methods
2. **ARI vs Depth** - Similar but using Adjusted Rand Index
3. **Per-cell-type stability** - Heatmap showing which cell types are stable at each depth
4. **Identify inflection point** - Where does agreement plateau?

---

## Directory Structure

```
downsampling-experiment/
├── main.nf                 # Nextflow pipeline
├── nextflow.config         # Configuration
├── bin/
│   ├── download_data.py    # Fetch data from CellxGene
│   ├── downsample.py       # Multinomial downsampling
│   ├── predict_seurat.R    # Seurat label transfer (from eval-refs)
│   ├── predict_scvi.py     # scVI label transfer (from eval-refs)
│   ├── compare_replicates.py  # Calculate agreement metrics
│   └── plot_stability.py   # Generate final plots
├── data/                   # Downloaded/processed data
├── results/                # Output directory
└── work/                   # Nextflow work directory
```

---

## Key Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| n_cells | 1000 | Number of cells to subsample |
| depth_levels | 100, 500, 1K, 2.5K, 5K, 10K, 25K, 50K, 100K | Reads per cell |
| n_replicates | 2 | Paired datasets at each depth |
| ref_key | "cell_type" | CellxGene cell type annotation field |
| seed | 42, 123 | Different seeds for each replicate |

---

## Expected Output

1. **Stability curve** showing agreement vs depth
2. **Minimum depth recommendation** for stable cell typing
3. **Method comparison** - does Seurat or scVI stabilize at different depths?
4. **Cell-type-specific thresholds** - some types may need more depth than others

---

## Questions to Resolve

1. **Reference dataset**: Should we use the full non-downsampled data as reference, or also downsample the reference? (Recommendation: keep reference at full depth)

2. **Smart-seq depth characteristics**: What is the typical read depth in the reference data? This determines our maximum downsampling range.

3. **Normalization**: Should downsampling happen before or after normalization? (Must be before - we downsample raw counts)
