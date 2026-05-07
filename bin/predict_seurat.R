#!/usr/bin/env Rscript
#
# Seurat-based label transfer for downsampling experiment.
# Transfers cell type labels from reference to downsampled query.
#

source("/space/grp/rschwartz/rschwartz/downsampling-experiment/bin/seurat_functions.R")
library(Seurat)
library(dplyr)
library(argparse)
set.seed(123)
options(future.globals.maxSize = 10 * 1024^3)  # 5 GB
parser <- ArgumentParser(description = "Seurat label transfer for downsampling experiment")
parser$add_argument("--ref_path", type="character", default="/space/grp/rschwartz/rschwartz/downsampling-experiment/results/seurat_processed/depth_original/seed_original/reference.rds",
                    help="Path to reference .rds file")
parser$add_argument("--query_path", type="character", default="/space/grp/rschwartz/rschwartz/downsampling-experiment/results/seurat_processed/depth_50000/seed_123/downsampled_50000_123.rds",
                    help="Path to downsampled query .rds file")
parser$add_argument("--output_dir", type="character", default=".",
                    help="Output directory")
parser$add_argument("--label_key", type="character", default="cell_type",
                    help="Column name for cell type labels")
parser$add_argument("--integration_method", type="character", default="pcaproject",
                    help="Integration method")
parser$add_argument("--dims", type="integer", default=50,
                    help="Number of PCA dimensions")
parser$add_argument("--max.features", type="integer", default=200,
                    help="Maximum number of features for anchors")
parser$add_argument("--k.anchor", type="integer", default=10,
                    help="Number of anchors")
parser$add_argument("--k.score", type="integer", default=30,
                    help="k for anchor scoring")
parser$add_argument("--k.weight", type="integer", default=50,
                    help="k for weight transfer")
parser$add_argument("--normalization_method", type="character", default="SCT",
                    help="Normalization method (SCT or LogNormalize)")

args <- parser$parse_args()
# Override for reproducibility
# Load reference
message("Loading reference: ", args$ref_path)
ref <- readRDS(args$ref_path)
query <- readRDS(args$query_path)
# Preprocess query with same method as reference
message("Preprocessing query...")
if (args$normalization_method == "SCT") {
    query <- SCTransform(query, verbose = FALSE)
} else {
    query <- NormalizeData(query, verbose = FALSE)
    query <- FindVariableFeatures(query, nfeatures = 2000, verbose = FALSE)
    query <- ScaleData(query, verbose = FALSE)
}
query <- RunPCA(query, npcs = args$dims, verbose = FALSE)

# Transfer labels
message("Finding transfer anchors...")
anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    normalization.method = args$normalization_method,
    npcs = args$dims,
    dims = 1:args$dims,
    reduction = args$integration_method,
    max.features = args$max.features,
    k.anchor = args$k.anchor,
    k.score = args$k.score,
    reference.reduction = "pca",
    recompute.residuals = FALSE
)

# Adjust k.weight if needed
k.weight <- min(args$k.weight, floor(nrow(anchors@anchors) / args$k.score))

message("Transferring labels...")
predictions <- tryCatch({
    TransferData(
        anchorset = anchors,
        refdata = ref@meta.data[[args$label_key]],
        weight.reduction = args$integration_method,
        k.weight = k.weight
    )
}, error = function(e) {
    message("Error in TransferData: ", e$message)
    # Extract k.weight suggestion from error message
    match <- regmatches(e$message, regexec("less than ([0-9]+)", e$message))
    if (length(match[[1]]) > 1) {
        k.weight.new <- as.integer(match[[1]][2]) - 1
        message("Adjusting k.weight to: ", k.weight.new)
        TransferData(
            anchorset = anchors,
            refdata = ref@meta.data[[args$label_key]],
            weight.reduction = args$integration_method,
            k.weight = k.weight.new
        )
    } else {
        stop(e)
    }
})

# Format output
prediction_scores <- predictions %>%
    as.data.frame() %>%
    select(-any_of(c("predicted.id", "prediction.score.max"))) %>%
    rename_with(~ gsub("prediction.score.", "", .))

# Get predicted labels
predicted_labels <- predictions$predicted.id

# Create output data frame
output_df <- data.frame(
    cell_id = colnames(query),
    predicted_label = predicted_labels,
    max_score = predictions$prediction.score.max,
    stringsAsFactors = FALSE
)

# Add original labels if available
if (args$label_key %in% colnames(query@meta.data)) {
    output_df$original_label <- query@meta.data[[args$label_key]]
}

# Save outputs
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

query_name <- basename(args$query_path) %>% gsub(".h5ad", "", .)
ref_name <- basename(args$ref_path) %>% gsub(".rds", "", .)

# Save predictions
pred_path <- file.path(args$output_dir, paste0(query_name, "_seurat_predictions.tsv"))
write.table(output_df, file = pred_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved predictions to: ", pred_path)

# Save prediction scores
scores_path <- file.path(args$output_dir, paste0(query_name, "_seurat_scores.tsv"))
write.table(prediction_scores, file = scores_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved scores to: ", scores_path)

message("Done!")
