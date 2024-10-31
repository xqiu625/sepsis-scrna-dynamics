```R
# module_score_analysis.R
# Calculate and visualize gene module scores across cell types
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)

# Define gene modules
GENE_MODULES <- list(
  "OXPHOS" = c(
    "NDUFA10", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7",
    "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFB10", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4",
    # ... [rest of OXPHOS genes]
  ),
  
  "Glycolysis" = c(
    "ENO1", "ENO2", "ENO3", "ENTPD5", "GAPDHS", "GAPDH", "GCK", "GPD1", "GPI",
    # ... [rest of Glycolysis genes]
  ),
  
  "MHC_ClassII" = c(
    "CD74", "HLA-DRA", "HLA-DQA1", "HLA-DQA2", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1",
    # ... [rest of MHC Class II genes]
  )
)

# Cell types to analyze
CELL_TYPES <- c(
  'CD14+ Mono', 'FCGR3A+ Mono', 'Platelet', 'B', 
  'CD4+ T', 'CD8+ T', 'NK', 'Erythroid precursor', 'DC'
)

#' Calculate module scores
#' @param seurat_obj Seurat object
#' @param module_name Name of the gene module
#' @param gene_list List of genes in the module
#' @return Updated Seurat object
calculate_module_score <- function(seurat_obj, module_name, gene_list) {
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(gene_list),
    name = module_name,
    search = TRUE
  )
  return(seurat_obj)
}

#' Create boxplot for a cell type and module
#' @param data Module score data
#' @param cell_type Cell type to plot
#' @param module_name Name of the module
#' @return ggplot object
create_module_boxplot <- function(data, cell_type, module_name) {
  # Filter data for cell type
  plot_data <- data %>%
    filter(CellType == cell_type) %>%
    select(condition, matches(module_name))
  
  # Perform statistical test
  stat.test <- plot_data %>%
    wilcox_test(as.formula(paste0(module_name, "~ condition"))) %>%
    adjust_pvalue() %>%
    add_significance("p.adj") %>%
    add_xy_position(fun = "max", x = "condition")
  
  # Create plot
  ggboxplot(
    plot_data,
    x = "condition",
    y = module_name,
    fill = "condition",
    xlab = cell_type,
    palette = c("#abdda4", "#d7191c", "#2b83ba")  # HC, NS, S colors
  ) +
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks.x = element_blank()
    )
}

#' Process module scores for all cell types
#' @param seurat_obj Seurat object
#' @param module_name Module name
#' @param output_dir Output directory
process_module <- function(seurat_obj, module_name, output_dir) {
  # Get visualization data
  df_viz <- seurat_obj[[c(paste0(module_name, "1"), "condition", "CellType")]]
  
  # Create plots for each cell type
  plots <- map(CELL_TYPES, ~create_module_boxplot(df_viz, .x, paste0(module_name, "1")))
  
  # Combine plots
  combined_plot <- plot_grid(plotlist = plots, ncol = length(CELL_TYPES))
  
  # Save plot
  ggsave(
    file.path(output_dir, paste0(module_name, "_scores.png")),
    combined_plot,
    width = 24,
    height = 4,
    dpi = 300
  )
  
  return(combined_plot)
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
  
  # Process each module
  results <- map(names(GENE_MODULES), function(module_name) {
    message("Processing module: ", module_name)
    
    # Calculate scores
    seurat_obj <- calculate_module_score(
      seurat_obj,
      module_name,
      GENE_MODULES[[module_name]]
    )
    
    # Create visualizations
    process_module(seurat_obj, module_name, output_dir)
  })
  
  # Save processed object
  saveRDS(seurat_obj, file.path(output_dir, "seurat_with_modules.rds"))
  
  return(results)
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/module_scores"
# )
```
