# umap_conditions.R
# Generate UMAP visualizations with exact condition colors
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)

# Define color schemes
CELL_COLORS <- c(
  "B" = '#0868ac',                    # blue
  "CD14+ Mono" = '#6a51a3',          # dark purple
  "CD4+ T" = '#33a02c',              # light green
  "CD8+ T" = '#969696',              # light grey
  "CMP" = '#6a51a3',                 # purple
  "DC" = '#006d2c',                  # dark green
  "Erythroid precursor" = '#ef3b2c', # red
  "FCGR3A+ Mono" = '#9e9ac8',        # light purple
  "Neutrophil" = '#fb9a99',          # pink
  "NK" = '#525252',                  # dark grey
  "Platelet" = '#ff7f00'             # orange
)

# Updated condition colors to exact specifications
CONDITION_COLORS <- c(
  "HC" = "#abdda4",    # specified green
  "NS_LS" = "#d7191c", # specified red
  "NS_ES" = "#fdae61", # specified orange
  "S" = "#2b83ba"      # specified blue
)

#' Add condition and time metadata to Seurat object
#' @param seurat_obj Seurat object
#' @return Seurat object with updated metadata
add_metadata <- function(seurat_obj) {
  # Get metadata
  metadata <- seurat_obj@meta.data
  
  # Add condition labels
  metadata$condition <- case_when(
    str_detect(metadata$cells, "^HC") ~ "HC",
    str_detect(metadata$cells, "^NS1") ~ "NS_LS",  # Late stage
    str_detect(metadata$cells, "^NS2") ~ "NS_ES",  # Early stage
    str_detect(metadata$cells, "^S") ~ "S"
  )
  
  # Update Seurat object
  seurat_obj@meta.data <- metadata
  
  return(seurat_obj)
}

#' Create UMAP plot
#' @param seurat_obj Seurat object
#' @param show_legend Whether to show legend
#' @return ggplot object
create_umap <- function(seurat_obj, show_legend = TRUE) {
  plot <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "condition",
    cols = CONDITION_COLORS,
    pt.size = 0.1
  )
  
  if (!show_legend) {
    plot <- plot + NoLegend()
  }
  
  return(plot)
}

#' Save plot with consistent parameters
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width
#' @param height Plot height
save_plot <- function(plot, filename, width = 6, height = 6) {
  png(
    filename,
    width = width * 300,
    height = height * 300,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  print(plot)
  dev.off()
}

#' Main execution function
#' @param input_path Path to Seurat object
#' @param output_dir Output directory
main <- function(input_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  seurat_obj <- readRDS(input_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Add metadata
  seurat_obj <- add_metadata(seurat_obj)
  
  # Create UMAP with legend
  umap_with_legend <- create_umap(seurat_obj, show_legend = TRUE)
  save_plot(
    umap_with_legend,
    file.path(output_dir, "umap_conditions_with_legend.png")
  )
  
  # Create UMAP without legend
  umap_no_legend <- create_umap(seurat_obj, show_legend = FALSE)
  save_plot(
    umap_no_legend,
    file.path(output_dir, "umap_conditions_no_legend.png")
  )
  
  # Save processed object
  saveRDS(seurat_obj, file.path(output_dir, "processed_seurat.rds"))
  
  return(seurat_obj)
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/umap"
# )
