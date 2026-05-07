
library(Seurat)
library(reticulate)
library(SeuratDisk)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
if (length(script_dir) == 0) script_dir <- "."
source(file.path(script_dir, "seurat_functions.R"))
options(future.globals.maxSize = 10 * 1024^3)  # 5 GB

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", default = "/space/grp/rschwartz/rschwartz/eval-references/work/9f/8e0b151b44c5e5f9c8bc54915c1040/whole_cortex.h5ad")
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="LogNormalize")
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=30)
parser$add_argument("--batch_key", type="character", help="Batch key", default="dataset_id")
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)
parser$add_argument("--batch_correct", type="character", help="Batch correct (TRUE/FALSE)", default="TRUE")
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=5)
parser$add_argument("--k.weight", type="integer", help="k.weight", default=30)
parser$add_argument("--k.score", type="integer", help="k.score", default=30)
args = parser$parse_args()
normalization_method = args$normalization_method
dims = args$dims
batch_key = args$batch_key
n_features = args$nfeatures
h5ad_file = args$h5ad_file
integration_method = args$integration_method
batch_correct = as.logical(args$batch_correct)
sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
k.anchor = args$k.anchor
k.weight = args$k.weight
k.score = args$k.score

batch_key_value <- eval(parse(text = paste0("sceasy_seurat$", batch_key)))

if ("feature_id" %in% colnames(sceasy_seurat@assays$RNA[[]])) {
  sceasy_seurat <- rename_features(sceasy_seurat, column_name="feature_id")
}

# If batch correction is enabled
if (batch_correct) {
    # Split by batch
    sceasy_list <- SplitObject(sceasy_seurat, split.by = batch_key)
    
    if (length(sceasy_list) == 1) {
        message("Only one batch found. Skipping integration.")
        batch_correct <- FALSE
    } else {
        # Filter out datasets with too few cells
        sceasy_list <- Filter(function(x) ncol(x) > dims, sceasy_list)

        if (normalization_method == "SCT") {
            assay <- "SCT"
            sceasy_list <- lapply(sceasy_list, function(x) SCTransform(x, verbose = FALSE, variable.features.n = n_features))
            features <- SelectIntegrationFeatures(sceasy_list)
            sceasy_list <- PrepSCTIntegration(sceasy_list, anchor.features = features)
        } else if (normalization_method == "LogNormalize") {
            assay <- "RNA"
            sceasy_list <- lapply(sceasy_list, function(x) {
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, nfeatures = n_features)
                x <- ScaleData(x)
                return(x)
            })
            features <- SelectIntegrationFeatures(sceasy_list)
        }

        # Find integration anchors
        anchors <- FindIntegrationAnchors(
            object.list = sceasy_list, dims = 1:dims, k.anchor = k.anchor, k.score = k.score,
            normalization.method = normalization_method, anchor.features = features
        )

        anchors_per_pair <- table(anchors@anchors$dataset1, anchors@anchors$dataset2)
        min_anchors <- min(anchors_per_pair[anchors_per_pair > 0])  # Ignore zero values
        k.anchor <- min(min_anchors, k.anchor)

        # Integrate data
        sceasy_seurat <- IntegrateData(
            anchorset = anchors, dims = 1:dims, k.weight = k.weight,
            normalization.method = normalization_method, new.assay.name = "integrated"
        )
        assay <- "integrated"
        if (normalization_method == "LogNormalize") {
            sceasy_seurat <- sceasy_seurat %>% ScaleData(assay=assay)
        }
    }
}

# If batch correction was skipped, perform normalization
if (!batch_correct) {
    if (normalization_method == "SCT") {
        sceasy_seurat <- SCTransform(sceasy_seurat, verbose = FALSE, variable.features.n = n_features)
        assay <- "SCT"
    } else if (normalization_method == "LogNormalize") {
        sceasy_seurat <- NormalizeData(sceasy_seurat, normalization.method = normalization_method)
        sceasy_seurat <- FindVariableFeatures(sceasy_seurat, nfeatures = n_features)
        sceasy_seurat <- ScaleData(sceasy_seurat)
        assay <- "RNA"
    }
}

# Run PCA on the appropriate assay
sceasy_seurat <- RunPCA(sceasy_seurat, npcs = dims, assay = assay)

filepath_prefix <- gsub(".h5ad", "", basename(h5ad_file))
saveRDS(sceasy_seurat, file = paste0(filepath_prefix, ".rds"))
new_prefix <- paste0(filepath_prefix, "_sct_integrated")
rds_to_h5ad(sceasy_seurat, filepath_prefix = new_prefix)
message(paste("Converted H5AD to RDS and saved to", paste0(filepath_prefix, ".rds")))
message(paste("Also saved integrated data to", paste0(new_prefix, ".h5ad")))