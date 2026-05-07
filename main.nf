#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Log parameters
log.info """
=========================================
Downsampling Experiment Pipeline
=========================================
Input data      : ${params.input_data}
Reference       : ${params.reference}
Output dir      : ${params.outdir}
N cells         : ${params.n_cells ?: 'all'}
Depth levels    : ${params.depth_levels}
Label key       : ${params.label_key}
Seeds           : ${params.seeds}
=========================================
"""

process DOWNLOAD_DATA {
    conda params.scanpy_env
    publishDir "${params.outdir}/downloaded", mode: 'copy'

    output:
    path "${params.output_prefix}.h5ad", emit: adata

    script:
    """
    python ${projectDir}/bin/download_data.py \
        --dataset_name "${params.dataset_name}" \
        --organism ${params.organism} \
        --census_version ${params.census_version} \
        ${params.n_cells != null ? "--n_cells ${params.n_cells} \\" : ""} \
        --output ${params.output_prefix}.h5ad
    """
}

process SPLIT_DATA {
    conda params.scanpy_env
    publishDir "${params.outdir}/data", mode: 'copy'

    input:
    path adata

    output:
    path "reference.h5ad", emit: reference
    path "query.h5ad", emit: query

    script:
    """
    python ${projectDir}/bin/split_data.py \\
        --input ${adata} \\
        --ref_output reference.h5ad \\
        --query_output query.h5ad \\
        --test_size ${params.test_size} \\
        --stratify_key ${params.stratify_key} \\
        --seed 42
    """
}


process DOWNSAMPLE {
    conda params.scanpy_env
    publishDir "${params.outdir}/downsampled/depth_${depth}/rep_${seed}", mode: 'copy'

    input:
    tuple val(depth), val(seed), path(adata)


    output:
    tuple val(depth), val(seed), path("downsampled_${depth}_${seed}.h5ad"), emit: downsampled
    path "downsampled_${depth}_${seed}_stats.tsv", emit: stats

    script:
    """
    python ${projectDir}/bin/downsample.py \\
        --input ${adata} \\
        --output downsampled_${depth}_${seed}.h5ad \\
        --target_depth ${depth} \\
        --seed ${seed}
    """
}

process SEURAT_PROCESSING {
    conda params.r_env
    publishDir "${params.outdir}/seurat_processed/depth_${depth}/seed_${seed}", mode: 'copy'

    input:
    tuple val(depth), val(seed), path(adata)

    output:
    tuple val(depth), val(seed), path("*.rds"), emit: seurat_rds

    script:
    """
    Rscript ${projectDir}/bin/seurat_processing.R \\
        --h5ad_file ${adata} \\
        --normalization_method ${params.normalization_method} \\
        --dims ${params.dims} \\
        --batch_correct FALSE
    """
}

process SEURAT_REF_PROCESSING {
    conda params.r_env
    publishDir "${params.outdir}/seurat_processed/reference", mode: 'copy'

    input:
    path adata

    output:
    path "*.rds", emit: seurat_rds

    script:
    """
    Rscript ${projectDir}/bin/seurat_processing.R \\
        --h5ad_file ${adata} \\
        --normalization_method ${params.normalization_method} \\
        --dims ${params.dims} \\
        --batch_correct FALSE
    """
}

process PREDICT_SEURAT {
    conda params.r_env
    publishDir "${params.outdir}/predictions/seurat/depth_${depth}/rep_${seed}", mode: 'copy'

    input:
    tuple val(depth), val(seed), path(query_rds), path(ref_rds)

    output:
    tuple val(depth), val(seed), val("seurat"), path("*_seurat_predictions.tsv"), emit: predictions

    script:
    """
    Rscript ${projectDir}/bin/predict_seurat.R \\
        --ref_path ${ref_rds} \\
        --query_path ${query_rds} \\
        --output_dir . \\
        --label_key ${params.label_key} \\
        --dims ${params.dims} \\
        --normalization_method ${params.normalization_method}
    """
}

process SETUP_MODEL {
    conda params.scanpy_env
    publishDir "${params.outdir}/model", mode: 'copy'

    output:
    path "scvi-*", emit: model_dir

    script:
    """
    python ${projectDir}/bin/setup.py \
        --organism ${params.organism} \
        --census_version ${params.census_version}
    """
}

process PREDICT_SCVI {
    conda params.scanpy_env
    publishDir "${params.outdir}/predictions/scvi/depth_${depth}/rep_${seed}", mode: 'copy'

    input:
    tuple val(depth), val(seed), path(query_h5ad), path(ref_h5ad), path(model_dir)

    output:
    tuple val(depth), val(seed), val("scvi"), path("*_scvi_predictions.tsv"), emit: predictions

    script:
    """
    python ${projectDir}/bin/predict_scvi.py \\
        --ref_path ${ref_h5ad} \\
        --query_path ${query_h5ad} \\
        --model_path ${model_dir} \\
        --output_dir . \\
        --label_key ${params.label_key} \\
        --seed ${seed}
    """
}

process COMBINE_DOWNSAMPLE_STATS {
    conda params.scanpy_env
    publishDir "${params.outdir}/downsampling_stats", mode: 'copy'

    input:
    path stats_files

    output:
    path "downsampling_stats.tsv"

    script:
    """
    # Combine all stats files into one, keeping header from first file only
    head -1 ${stats_files[0]} > downsampling_stats.tsv
    for f in ${stats_files}; do
        tail -n +2 "\$f" >> downsampling_stats.tsv
    done
    """
}

process COMPARE_REPLICATES {
    conda params.scanpy_env
    publishDir "${params.outdir}/comparisons/${method}/depth_${depth}", mode: 'copy'

    input:
    tuple val(depth), val(method), val(seed1), val(seed2), path(rep1_pred), path(rep2_pred)

    output:
    path "comparison_${method}_${depth}.tsv", emit: metrics
    path "comparison_${method}_${depth}_per_type.tsv", emit: per_type

    script:
    """
    python ${projectDir}/bin/compare_replicates.py \\
        --rep1 ${rep1_pred} \\
        --rep2 ${rep2_pred} \\
        --output comparison_${method}_${depth}.tsv \\
        --depth ${depth} \\
        --method ${method} \\
        --seed1 ${seed1} \\
        --seed2 ${seed2}
    """
}

process PLOT_RESULTS {
    conda params.scanpy_env
    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    path comparison_files

    output:
    path "*.png"
    path "*.tsv"

    script:
    """
    # Create input directory with all comparison files
    mkdir -p input_metrics
    cp ${comparison_files} input_metrics/

    python ${projectDir}/bin/plot_stability.py \\
        --input_dir input_metrics \\
        --output_dir . \\
        --stability_threshold ${params.stability_threshold}
    """
}

workflow {
    // Step 1: Download or use existing data
    if (params.input_data) {
        data_ch = Channel.fromPath(params.input_data)
    } else {
        DOWNLOAD_DATA()
        data_ch = DOWNLOAD_DATA.out.adata
    }

    // Step 2: Split data into reference and query
    SPLIT_DATA(data_ch)
    ref_ch = SPLIT_DATA.out.reference
    query_ch = SPLIT_DATA.out.query

    // Step 3: Download pre-trained scVI model
    SETUP_MODEL()
    model_ch = SETUP_MODEL.out.model_dir

    // Step 4: Create downsampled versions at each depth with different seeds
    depth_ch = Channel.from(params.depth_levels)
    seed_ch = Channel.from(params.seeds)

    // Combine depths and seeds with the query data
    downsample_input = depth_ch
        .combine(seed_ch)
        .combine(query_ch)

    DOWNSAMPLE(downsample_input)
    downsampled_ch = DOWNSAMPLE.out.downsampled

    // Combine all downsampling stats into one file
    COMBINE_DOWNSAMPLE_STATS(DOWNSAMPLE.out.stats.collect())

    // Step 5a: Process reference data with Seurat
    SEURAT_REF_PROCESSING(ref_ch)
    ref_rds_ch = SEURAT_REF_PROCESSING.out.seurat_rds

    // Step 5b: Process downsampled query data with Seurat
    SEURAT_PROCESSING(downsampled_ch)
    query_rds_ch = SEURAT_PROCESSING.out.seurat_rds

    // Seurat predictions: combine query RDS with reference RDS
    seurat_input = query_rds_ch
        .combine(ref_rds_ch)
        .map { depth, seed, query_rds, ref_rds ->
            tuple(depth, seed, query_rds, ref_rds)
        }

    PREDICT_SEURAT(seurat_input)

    // scVI predictions: same tuple structure (depth, seed, query_h5ad, ref_h5ad)
    scvi_input = downsampled_ch
        .combine(ref_ch)
        .combine(model_ch)
        .map { depth, seed, query_h5ad, ref_h5ad, model_dir ->
            tuple(depth, seed, query_h5ad, ref_h5ad, model_dir)
        }

    PREDICT_SCVI(scvi_input)

    // Step 6: Group replicates by depth and method for comparison
    all_predictions = PREDICT_SCVI.out.predictions
        .mix(PREDICT_SEURAT.out.predictions)

    // Group by depth and method, expecting 2 replicates each
    grouped_predictions = all_predictions
        .map { depth, seed, method, pred_file ->
            [depth, method, seed, pred_file]
        }
        .groupTuple(by: [0, 1], size: 2)
        .map { depth, method, seeds, pred_files ->
            [depth, method, seeds[0], seeds[1], pred_files[0], pred_files[1]]
        }

    COMPARE_REPLICATES(grouped_predictions)

    // Step 6: Collect all comparison results and plot
    all_comparisons = COMPARE_REPLICATES.out.metrics
        .mix(COMPARE_REPLICATES.out.per_type)
        .collect()

    PLOT_RESULTS(all_comparisons)
}

workflow.onComplete {
    log.info """
=========================================
Pipeline Complete!
=========================================
Duration    : ${workflow.duration}
Success     : ${workflow.success}
Results dir : ${params.outdir}/results
=========================================
"""
}

workflow.onError {
    log.error "Pipeline failed"
}
