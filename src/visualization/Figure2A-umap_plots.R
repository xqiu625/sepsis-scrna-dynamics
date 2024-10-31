# umap_visualization.R
# Generate UMAP visualizations for sepsis single-cell data
# Data from Qiu et al., 2021

library(Seurat)

#' Define cell type mapping for cluster annotation
#' @return Named vector of cluster to cell type mappings
get_cell_type_mapping <- function() {
  c(
    "0" = "B",
    "1" = "CD14+ Mono",
    "2" = "CD4+ T",
    "3" = "FCGR3A+ Mono",
    "4" = "CD4+ T",
    "5" = "CD14+ Mono",
    "6" = "CD4+ T",
    "7" = "CD8+ T",
    "8" = "NK",
    "9" = "B",
    "10" = "CD4+ T",
    "11" = "NK",
    "12" = "FCGR3A+ Mono",
    "13" = "DC",
    "14" = "B",
    "15" = "Platelet",
    "16" = "Platelet",
    "17" = "CD14+ Mono",
    "18" = "B",
    "19" = "B",
    "20" = "CD8+ T",
    "21" = "Erythroblast",
    "22" = "FCGR3A+ Mono",
    "23" = "CD4+ T",
    "24" = "B",
    "25" = "Platelet",
    "26" = "Neutrophil",
    "27" = "CD14+ Mono",
    "28" = "CD14+ Mono",
    "29" = "CD8+ T",
    "30" = "B",
    "31" = "CD4+ T",
    "32" = "NK",
    "33" = "CMP",
    "34" = "NK"
  )
}

#' Create UMAP visualization
#' @param seurat_obj Seurat object
#' @param split_by Variable to split by (optional)
#' @param show_label Whether to show labels
#' @param show_legend Whether to show legend
#' @param ncol Number of columns for split plots
#' @return UMAP plot
create_umap <- function(seurat_obj, 
                       split_by = NULL,
                       show_label = FALSE,
                       show_legend = TRUE,
                       ncol = NULL) {
  
  plot <- DimPlot(
    seurat_obj,
    reduction = "umap",
    split.by = split_by,
    label = show_label,
    label.size = 4,
    pt.size = 0.5,
    ncol = ncol
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
save_plot <- function(plot, filename, width = 8, height = 6) {
  dpi <- 300
  png(
    file = filename,
    width = dpi * width,
    height = dpi * height,
    units = "px",
    res = dpi,
    type = 'cairo'
  )
  print(plot)
  dev.off()
}

#' Generate all UMAP visualizations
#' @param input_path Path to Seurat object
#' @param output_dir Output directory
main <- function(input_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  seurat_obj <- readRDS(input_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Rename cluster identities
  seurat_obj <- RenameIdents(
    object = seurat_obj,
    get_cell_type_mapping()
  )
  
  # Generate plots
  
  # 1. Basic UMAP
  basic_umap <- create_umap(seurat_obj)
  save_plot(
    basic_umap,
    file.path(output_dir, "umap_celltype.png"),
    width = 8,
    height = 6
  )
  
  # 2. UMAP split by sample with labels
  split_umap_labeled <- create_umap(
    seurat_obj,
    split_by = "sample",
    show_label = TRUE,
    show_legend = FALSE,
    ncol = 4
  )
  save_plot(
    split_umap_labeled,
    file.path(output_dir, "umap_celltype_sample_labeled.png"),
    width = 24,
    height = 24
  )
  
  # 3. UMAP split by sample without labels
  split_umap <- create_umap(
    seurat_obj,
    split_by = "sample",
    ncol = 4
  )
  save_plot(
    split_umap,
    file.path(output_dir, "umap_celltype_sample.png"),
    width = 24,
    height = 24
  )
  
  # Save annotated object
  saveRDS(seurat_obj, 
          file = file.path(output_dir, "labeled_seurat.rds"))
  
  message("Visualizations complete. Results saved in: ", output_dir)
}

# Example usage:
# main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/visualizations"
# )
