# hif1a_heatmap.R
# Create heatmap of average HIF1A expression across conditions and cell types
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Define cell types in desired order
CELL_TYPES <- c(
  "B", "CD14+ Mono", "CD4+ T", "CD8+ T", "CMP", "DC",
  "Erythroid precursor", "FCGR3A+ Mono", "Neutrophil", "NK", "Platelet"
)

#' Calculate average gene expression for a condition
#' @param seurat_obj Seurat object
#' @param gene Gene name
#' @return Data frame of average expression
calculate_average_expression <- function(seurat_obj, gene = "HIF1A") {
  avg_expr <- AverageExpression(
    seurat_obj,
    group.by = "CellType",
    slot = "data",
    return.seurat = TRUE
  )
  
  # Extract and process expression data
  expr_df <- avg_expr@assays$RNA@data %>%
    as.data.frame() %>%
    subset(rownames(.) == gene) %>%
    select(all_of(CELL_TYPES))
  
  return(expr_df)
}

#' Normalize matrix to [-1,1] range
#' @param x Numeric matrix
#' @return Normalized matrix
range_normalize <- function(x) {
  2 * ((x - min(x)) / (max(x) - min(x))) - 1
}

#' Create expression heatmap
#' @param expr_matrix Expression matrix
#' @param output_file Output file path
create_heatmap <- function(expr_matrix, output_file) {
  # Scale and normalize data
  scaled_matrix <- scale(expr_matrix)
  normalized_matrix <- range_normalize(scaled_matrix)
  
  # Create color palette
  colors <- colorRampPalette(
    rev(brewer.pal(n = 7, name = "Spectral"))
  )(200)
  
  # Create and save heatmap
  png(
    output_file,
    width = 4 * 300,
    height = 3.2 * 300,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  
  pheatmap(
    normalized_matrix,
    border_color = NA,
    color = colors,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 10,
    main = "HIF1A Expression"
  )
  
  dev.off()
}

#' Process all conditions and create heatmap
#' @param seurat_obj Seurat object
#' @param output_dir Output directory
main <- function(seurat_obj, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define conditions
  conditions <- c("HC", "NS LS", "NS ES", "S")
  
  # Calculate expression for each condition
  expression_data <- map_df(conditions, function(cond) {
    subset_obj <- subset(seurat_obj, condition == cond)
    calculate_average_expression(subset_obj)
  }, .id = "condition")
  
  # Format matrix for heatmap
  heatmap_matrix <- expression_data %>%
    as.data.frame() %>%
    column_to_rownames("condition")
  
  # Create heatmap
  create_heatmap(
    heatmap_matrix,
    file.path(output_dir, "hif1a_expression_heatmap.png")
  )
  
  # Save processed data
  write.csv(
    heatmap_matrix,
    file.path(output_dir, "hif1a_expression_data.csv")
  )
  
  return(heatmap_matrix)
}

# Example usage:
# results <- main(
#   seurat_obj = seurat_subset,
#   output_dir = "results/hif1a_analysis"
# )
