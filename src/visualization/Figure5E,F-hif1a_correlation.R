# hif1a_correlation.R
# Analyze correlations between HIF1A and metabolic pathways
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)

# Define gene sets
GENE_SETS <- list(
  OXPHOS = c(
    "ACTN3", "AK2", "ATP5F1C", "ATP5F1D", "ATP7A", "CHCHD10", 
    "COQ7", "COQ9", "COX10", "COX15", "COX4I1", "COX5A", "COX8A",
    # ... [rest of OXPHOS genes]
  ),
  
  Glycolysis = c(
    "ALDOA", "ALDOC", "ARNT", "ENO1", "ENO2", "ENO3", "ENTPD5",
    "GAPDH", "GAPDHS", "GCK", "GPI", "HIF1A", "HK1", "HK2", "HK3",
    # ... [rest of Glycolysis genes]
  )
)

# Cell types to analyze
CELL_TYPES <- c(
  'CD14+ Mono', 'FCGR3A+ Mono', 'Platelet', 'B', 'CD4+ T',
  'CD8+ T', 'NK', 'Erythroid precursor', 'DC', 'Neutrophil'
)

# Color scheme for conditions
CONDITION_COLORS <- c(
  "HC" = "#4daf4a",
  "NS LS" = "#d95f0e",
  "NS ES" = "#e41a1c",
  "S" = "#377eb8"
)

#' Add condition metadata to Seurat object
#' @param seurat_obj Seurat object
#' @return Seurat object with condition metadata
add_condition_metadata <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data
  
  metadata$condition <- case_when(
    str_detect(metadata$cells, "^HC") ~ "HC",
    str_detect(metadata$cells, "^NS1") ~ "NS LS",
    str_detect(metadata$cells, "^NS2") ~ "NS ES",
    str_detect(metadata$cells, "^S") ~ "S"
  )
  
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}

#' Calculate module scores for gene sets
#' @param seurat_obj Seurat object
#' @return Seurat object with module scores
calculate_module_scores <- function(seurat_obj) {
  for (name in names(GENE_SETS)) {
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(GENE_SETS[[name]]),
      name = name,
      search = TRUE
    )
  }
  return(seurat_obj)
}

#' Create correlation plot
#' @param data Filtered data for one cell type
#' @param y_var Variable to correlate with HIF1A
#' @param cell_type Cell type name for title
#' @return ggplot object
create_correlation_plot <- function(data, y_var, cell_type) {
  ggscatter(
    data,
    x = "HIF1A",
    y = y_var,
    color = "Condition",
    size = 0.1,
    add = "reg.line",
    conf.int = TRUE,
    palette = CONDITION_COLORS
  ) +
    stat_cor(aes(color = Condition), label.x = 1.5) +
    labs(title = cell_type) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

#' Process correlations for all cell types
#' @param data Processed expression data
#' @param y_var Variable to correlate with HIF1A
#' @return List of plots
process_correlations <- function(data, y_var) {
  map(CELL_TYPES, function(cell_type) {
    data %>%
      filter(CellType == cell_type) %>%
      create_correlation_plot(y_var, cell_type)
  })
}

#' Save correlation plots
#' @param plots List of plots
#' @param filename Output filename
save_correlation_plots <- function(plots, filename) {
  png(
    filename,
    width = 24 * 300,
    height = 12 * 300,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  
  print(plot_grid(
    plotlist = plots,
    ncol = 5
  ))
  
  dev.off()
}

#' Main execution function
#' @param input_path Path to Seurat object
#' @param output_dir Output directory
main <- function(input_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load and process data
  seurat_obj <- readRDS(input_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Add metadata and calculate scores
  seurat_obj <- add_condition_metadata(seurat_obj)
  seurat_obj <- calculate_module_scores(seurat_obj)
  
  # Extract and process data
  expression_data <- FetchData(
    seurat_obj,
    c("OXPHOS1", "Glycolysis1", "condition", "CellType", "HIF1A")
  ) %>%
    as.data.frame() %>%
    filter(HIF1A != 0) %>%
    rename(
      OXPHOS = OXPHOS1,
      Glycolysis = Glycolysis1,
      Condition = condition
    )
  
  # Generate correlation plots
  oxphos_plots <- process_correlations(expression_data, "OXPHOS")
  glycolysis_plots <- process_correlations(expression_data, "Glycolysis")
  
  # Save plots
  save_correlation_plots(
    oxphos_plots,
    file.path(output_dir, "hif1a_oxphos_correlation.png")
  )
  
  save_correlation_plots(
    glycolysis_plots,
    file.path(output_dir, "hif1a_glycolysis_correlation.png")
  )
  
  return(list(
    data = expression_data,
    oxphos_plots = oxphos_plots,
    glycolysis_plots = glycolysis_plots
  ))
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/hif1a_correlation"
# )
