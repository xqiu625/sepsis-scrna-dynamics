# hladr_monocyte_analysis.R
# Analyze HLA-DR expression in monocyte subsets
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(ggpubr)

# Define constants
HLA_DR_GENES <- list(
  c("HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5")
)

MONOCYTE_TYPES <- c("CD14+ Mono", "FCGR3A+ Mono")

#' Calculate HLA-DR score for a cell subset
#' @param seurat_obj Seurat object
#' @param cell_type Cell type to analyze
#' @return Seurat object with HLA-DR score
calculate_hladr_score <- function(seurat_obj, cell_type) {
  # Subset data
  subset_obj <- subset(seurat_obj, CellType == cell_type)
  
  # Calculate score
  subset_obj <- AddModuleScore(
    object = subset_obj,
    features = HLA_DR_GENES,
    name = 'HLA-DR',
    search = TRUE
  )
  
  # Rename score column
  colnames(subset_obj@meta.data)[colnames(subset_obj@meta.data) == "HLA.DR1"] <- "HLA-DR"
  
  return(subset_obj)
}

#' Create violin plot for HLA-DR expression
#' @param seurat_obj Seurat object
#' @param output_file Output file path
create_hladr_plot <- function(seurat_obj, output_file) {
  # Create base plot
  plot <- VlnPlot(
    seurat_obj,
    features = "HLA-DR",
    sort = 'increasing',
    group.by = "sample_type",
    adjust = 1,
    pt.size = 0,
    combine = TRUE
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    ) +
    ggtitle("")
  
  # Save plot
  ggsave(
    output_file,
    plot = plot,
    width = 4,
    height = 4,
    dpi = 300
  )
  
  return(plot)
}

#' Main execution function
#' @param input_path Path to Seurat object
#' @param output_dir Output directory
main <- function(input_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  message("Loading Seurat object...")
  seurat_obj <- readRDS(input_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Add metadata
  message("Adding metadata...")
  seurat_obj <- add_detailed_metadata(seurat_obj)
  seurat_obj$CellType <- Idents(seurat_obj)
  
  # Process each monocyte type
  results <- map(MONOCYTE_TYPES, function(mono_type) {
    message("Processing ", mono_type, "...")
    
    # Calculate HLA-DR score
    mono_obj <- calculate_hladr_score(seurat_obj, mono_type)
    
    # Create plot
    output_file <- file.path(
      output_dir,
      paste0(gsub(" ", "_", mono_type), "_HLA-DR.png")
    )
    
    plot <- create_hladr_plot(mono_obj, output_file)
    
    return(list(
      data = mono_obj,
      plot = plot
    ))
  })
  
  names(results) <- MONOCYTE_TYPES
  return(results)
}

# Example usage:
# results <- main(
#   input_path = "data/sepsis_singlecell_transcriptome.rds",
#   output_dir = "results/hladr_analysis"
# )
