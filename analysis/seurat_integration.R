# seurat_integration.R
# Integration of sepsis single-cell data using Seurat v3
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)

preprocess_seurat <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000
  )
  return(seurat_obj)
}

integrate_samples <- function(seurat_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load merged object
  merged_seurat <- readRDS(seurat_path)
  
  # Split object by sample
  sample_list <- SplitObject(merged_seurat, split.by = "orig.ident")
  
  # Normalize and identify variable features for each dataset
  sample_list <- lapply(X = sample_list, FUN = preprocess_seurat)
  
  # Select features for integration
  features <- SelectIntegrationFeatures(
    object.list = sample_list,
    nfeatures = 2000
  )
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(
    object.list = sample_list,
    anchor.features = features,
    dims = 1:30
  )
  
  # Integrate data
  seurat_integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:30
  )
  
  DefaultAssay(seurat_integrated) <- "integrated"
  
  seurat_integrated <- ScaleData(seurat_integrated)
  seurat_integrated <- RunPCA(seurat_integrated, npcs = 50)
  
  pdf(file.path(output_dir, "elbow_plot.pdf"))
  print(ElbowPlot(seurat_integrated, ndims = 50))
  dev.off()
  
  # Run UMAP and clustering
  seurat_integrated <- RunUMAP(
    seurat_integrated,
    reduction = "pca",
    dims = 1:30
  )
  
  seurat_integrated <- FindNeighbors(
    seurat_integrated,
    reduction = "pca",
    dims = 1:30
  )
  
  seurat_integrated <- FindClusters(
    seurat_integrated,
    resolution = 1.0
  )
  
  # Save integrated object
  saveRDS(seurat_integrated, 
          file = file.path(output_dir, "labeled_clustered_integrated_sepsis.rds"))
  
  # Generate QC plots
  # UMAP by sample
  pdf(file.path(output_dir, "umap_by_sample.pdf"))
  print(DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident"))
  dev.off()
  
  # UMAP by condition
  pdf(file.path(output_dir, "umap_by_condition.pdf"))
  print(DimPlot(seurat_integrated, reduction = "umap", group.by = "condition"))
  dev.off()
  
  # UMAP by clusters
  pdf(file.path(output_dir, "umap_by_clusters.pdf"))
  print(DimPlot(seurat_integrated, reduction = "umap", label = TRUE))
  dev.off()
  
  return(seurat_integrated)
}

# Add condition and time labels
add_condition_labels <- function(seurat_obj) {
  # Add condition labels
  seurat_obj$condition <- NA
  seurat_obj$condition[which(str_detect(colnames(seurat_obj), "^HC"))] <- "HC"
  seurat_obj$condition[which(str_detect(colnames(seurat_obj), "^NS_LS"))] <- "NS_LS"
  seurat_obj$condition[which(str_detect(colnames(seurat_obj), "^NS_ES"))] <- "NS_ES"
  seurat_obj$condition[which(str_detect(colnames(seurat_obj), "^S"))] <- "S"
  
  # Add time point labels
  seurat_obj$time <- NA
  seurat_obj$time[which(str_detect(colnames(seurat_obj), "_T0"))] <- "T0"
  seurat_obj$time[which(str_detect(colnames(seurat_obj), "_T6"))] <- "T6"
  seurat_obj$time[which(str_detect(colnames(seurat_obj), "^HC"))] <- "NA"
  
  return(seurat_obj)
}

# Define canonical markers for cell type annotation
canonical_markers <- list(
  "B" = c("MS4A1"),
  "CD14+ Mono" = c("CD14", "LYZ"),
  "CD4+ T" = c("IL7R", "CCR7", "CD27"),
  "CD8+ T" = c("CD8A"),
  "DC" = c("FCER1A", "CST3", "CD123", "GZMB"),
  "Erythroid precursor" = c("GYPB", "AHSP"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
  "Neutrophil" = c("JAML", "SERPINB"),
  "NK" = c("GNLY", "NKG7"),
  "Platelet" = c("PPBP")
)

# Annotate cell types
annotate_cell_types <- function(seurat_obj, markers = canonical_markers) {
  # Find markers for each cluster
  DefaultAssay(seurat_obj) <- "RNA"
  all.markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  # Save marker results
  write.csv(all.markers, "results/cluster_markers.csv")
  
  # Generate expression plots for canonical markers
  marker_genes <- unlist(markers)
  pdf("results/canonical_marker_expression.pdf", width = 12, height = 8)
  print(FeaturePlot(seurat_obj, features = marker_genes, ncol = 4))
  dev.off()
  
  return(seurat_obj)
}

# Main execution
main <- function() {
  # Set paths
  input_path <- "results/qc/sepsis_merged_filtered.rds"
  output_dir <- "results/integrated"
  
  # Run integration
  integrated_seurat <- integrate_samples(input_path, output_dir)
  
  # Add metadata
  integrated_seurat <- add_condition_labels(integrated_seurat)
  
  # Annotate cell types
  integrated_seurat <- annotate_cell_types(integrated_seurat)
  
  # Save final object
  saveRDS(integrated_seurat, 
          file = "results/integrated/labeled_clustered_integrated_sepsis.rds")
}

# Run pipeline
main()
