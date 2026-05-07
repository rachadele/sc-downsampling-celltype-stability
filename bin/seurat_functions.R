load_required_packages <- function() {
  packages <- c("Seurat", "dplyr", "ggplot2", "data.table", "readxl", "stringr", "harmony","leiden", 
  "gprofiler2", "GEOquery", "ggExtra","SeuratDisk","tidyr","patchwork","RColorBrewer")

  for (pkg in packages) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
}

rename_features <- function(seurat_obj, column_name) {
    counts <- seurat_obj[["RNA"]]@counts
    data <- seurat_obj[["RNA"]]@data 

    feature_meta <- seurat_obj[["RNA"]][[]]
    feature_meta[["orig_features"]] <- rownames(feature_meta)
    rownames(feature_meta) <- feature_meta[[column_name]]

    rownames(counts) <- rownames(feature_meta)
    rownames(data) <- rownames(feature_meta)

    newRNA <- CreateAssayObject(counts=counts)
    newRNA$data <- data
    newRNA[[]] <- feature_meta

   seurat_obj[["RNA"]] <- newRNA
   DefaultAssay(seurat_obj) <- "RNA"
   return(seurat_obj)
}

rds_to_h5ad <- function(obj, filepath_prefix) {
    SeuratDisk::SaveH5Seurat(obj, filename = paste0(filepath_prefix,".h5Seurat"), overwrite = TRUE)
    SeuratDisk::Convert(paste0(filepath_prefix,".h5Seurat"), dest = "h5ad")
}

h5ad_to_rds <- function(h5ad_file_path){
    SeuratDisk::Convert(h5ad_file_path, dest = "h5seurat", overwrite = FALSE)
    message("Reading H5Seurat...")
    h5seurat_file_path <- gsub(".h5ad", ".h5seurat", h5ad_file_path)
    seurat_obj <- SeuratDisk::LoadH5Seurat(h5seurat_file_path, assays = "RNA")
    message("Read Seurat object:")
    return(seurat_obj)
}

transfer_multiple_labels <- function(
            query, reference, reduction, normalization_method="SCT",
            ref_keys, dims, max.features, 
            k.anchor, k.score, k.weight, project.query=NULL) {
    
    #threshold_dims <- 2000
    if (is.null(project.query)) {
        if ((ncol(query) > ncol(reference))) {
            project.query = TRUE
        } else {
            project.query = FALSE
        }
    }
    anchors <- FindTransferAnchors(reference=reference, 
        normalization.method=normalization_method, 
        query=query, 
        # assay should default to whatever normalization was applied
        npcs=dims, dims=1:dims, 
        reduction = reduction, # had to remove project query for SCTransform to work
        max.features=max.features, 
        k.anchor=k.anchor, 
        k.score=k.score,
        reference.reduction="pca", # use precomputed PCA from seurat_processing step
        recompute.residuals=FALSE) # setting to false fixes some errors

    k.weight = min(k.weight, floor(nrow(anchors@anchors) / k.score ))
    key = ref_keys[1] # assumes keys are ordered

    predictions <- tryCatch({
        TransferData(anchorset = anchors, refdata=key, reference=reference, weight.reduction=reduction, k.weight = k.weight)
        }, error = function(e) {
        
        message("Error in TransferData: ", e$message)
        # Extract k.weight suggestion (e.g., "less than 27") from error message
        match <- regmatches(e$message, regexec("less than ([0-9]+)", e$message))
        k.weight.new <- as.integer(match[[1]][2]) - 1 # Decrease by 1 to ensure it is less than the suggested value
        message("Adjusting k.weight to: ", k.weight.new)
        # Retry TransferData with the corrected k.weight
        TransferData(
            anchorset = anchors,
            refdata = key,
            reference = reference,
            weight.reduction = reduction,
            k.weight = k.weight.new
        )
    })   
    return(predictions)

}