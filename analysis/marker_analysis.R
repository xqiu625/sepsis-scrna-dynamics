# marker_analysis.R
# Differential expression analysis and cell type annotation
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(MAST)
library(scCATCH)

#' Perform differential expression analysis using MAST
#' @param seurat_obj Seurat object
#' @param output_dir Output directory
#' @return Seurat object with marker genes
find_marker_genes <- function(seurat_obj, output_dir) {
  # Set RNA as default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find markers using MAST
  markers <- FindAllMarkers(
    object = seurat_obj,
    assay = 'RNA',
    only.pos = TRUE,
    test.use = 'MAST'
  )
  
  # Store markers in object
  seurat_obj@misc$markers <- markers
  
  # Save marker results
  write.table(
    markers,
    file = file.path(output_dir, "marker_genes_mast.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(seurat_obj)
}

#' Generate marker gene heatmap
#' @param seurat_obj Seurat object
#' @param markers Marker gene dataframe
#' @param output_dir Output directory
plot_marker_heatmap <- function(seurat_obj, markers, output_dir) {
  # Get top 10 markers per cluster
  top10 <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = 10, wt = avg_logFC)
  
  # Generate heatmap
  heatmap <- DoHeatmap(
    object = subset(seurat_obj, downsample = 500),
    features = top10$gene
  ) + 
    NoLegend()
  
  # Save heatmap
  ggsave(
    file = file.path(output_dir, "marker_heatmap.pdf"),
    plot = heatmap,
    width = 20,
    height = 16,
    units = "in",
    limitsize = FALSE
  )
}

#' Run scCATCH analysis for cell type annotation
#' @param seurat_obj Seurat object
#' @return Seurat object with scCATCH results
run_sccatch <- function(seurat_obj) {
  # Find marker genes
  sccatch_markers <- findmarkergenes(
    seurat_obj,
    species = 'Human',
    cluster = 'All',
    match_CellMatch = TRUE,
    cancer = NULL,
    tissue = "Blood",
    cell_min_pct = 0.25,
    logfc = 0.25,
    pvalue = 0.05
  )
  
  # Annotate cell types
  sccatch_celltype <- scCATCH(
    object = sccatch_markers,
    species = 'Human',
    cancer = NULL,
    tissue = "Blood"
  )
  
  # Store results in object
  seurat_obj@misc$scCATCH_markers <- sccatch_markers
  seurat_obj@misc$scCATCH_cellType <- sccatch_celltype
  
  return(seurat_obj)
}

#' Main execution function
#' @param input_path Path to Seurat object
#' @param output_dir Output directory
main <- function(input_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  seurat_obj <- readRDS(input_path)
  
  # Find marker genes
  seurat_obj <- find_marker_genes(seurat_obj, output_dir)
  
  # Plot heatmap
  plot_marker_heatmap(seurat_obj, seurat_obj@misc$markers, output_dir)
  
  # Run scCATCH
  seurat_obj <- run_sccatch(seurat_obj)
  
  # Save annotated object
  saveRDS(seurat_obj, file = file.path(output_dir, "annotated_seurat.rds"))
  
  return(seurat_obj)
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/markers"
# )
