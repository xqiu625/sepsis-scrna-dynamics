# marker_analysis_scCATCH.R
# Marker gene identification and cell type annotation
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
  
  # Save marker results
  write.table(
    markers,
    file = file.path(output_dir, "marker_genes_mast.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(list(
    seurat_obj = seurat_obj,
    markers = markers
  ))
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

#' Perform cell type annotation using scCATCH
#' @param markers Marker gene dataframe
#' @param output_dir Output directory
#' @return Cell type annotations
annotate_cell_types <- function(markers, output_dir) {
  # Run scCATCH annotation
  cell_types <- scCATCH(
    object = markers,
    species = 'Human',
    cancer = NULL,
    tissue = 'Blood'
  )
  
  # Save annotations
  write.table(
    cell_types,
    file = file.path(output_dir, "cell_type_annotations.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(cell_types)
}

#' Generate annotation summary
#' @param annotations Cell type annotations
#' @param output_dir Output directory
create_annotation_summary <- function(annotations, output_dir) {
  # Create summary
  summary <- annotations %>%
    group_by(cell_type) %>%
    summarise(
      clusters = paste(cluster, collapse = ","),
      n_markers = n(),
      top_markers = paste(head(gene, 5), collapse = ",")
    )
  
  # Save summary
  write.table(
    summary,
    file = file.path(output_dir, "annotation_summary.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(summary)
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
  results <- find_marker_genes(seurat_obj, output_dir)
  
  # Plot heatmap
  plot_marker_heatmap(results$seurat_obj, results$markers, output_dir)
  
  # Annotate cell types
  annotations <- annotate_cell_types(results$markers, output_dir)
  
  # Create summary
  summary <- create_annotation_summary(annotations, output_dir)
  
  # Add annotations to Seurat object
  results$seurat_obj@misc$cell_types <- annotations
  
  # Save annotated object
  saveRDS(results$seurat_obj, 
          file = file.path(output_dir, "annotated_seurat.rds"))
  
  message("Analysis complete. Results saved in: ", output_dir)
  
  return(list(
    seurat_obj = results$seurat_obj,
    markers = results$markers,
    annotations = annotations,
    summary = summary
  ))
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/markers"
# )
