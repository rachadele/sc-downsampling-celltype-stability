#!/bin/bash
# Run the downsampling experiment pipeline

# Human Smart-seq data
nextflow run main.nf \
    --dataset_name "Multiple Cortical Areas SMART-seq" \
    --organism homo_sapiens \
    --label_key cell_type \
    -resume

# Mouse data (uncomment to use)
# nextflow run main.nf \
#     --dataset_name "your_mouse_dataset_name" \
#     --organism mus_musculus \
#     --n_cells 1000 \
#     --depth_levels '[10000, 50000, 100000, 250000, 500000, 1000000]' \
#     --label_key cell_type \
#     --outdir results_mouse \
#     -resume

# Use existing h5ad (uncomment to use)
# nextflow run main.nf \
#     --input_data /path/to/your.h5ad \
#     --outdir results \
#     -resume
