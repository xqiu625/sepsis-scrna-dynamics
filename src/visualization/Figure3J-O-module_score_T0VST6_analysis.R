# module_score_analysis.R
# Calculate and visualize gene module scores across conditions and cell types
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(rstatix)

# Define gene modules
GENE_MODULES <- list(
  OXPHOS = c("NDUFA10", "NDUFA1", "NDUFA2", /* ... */),
  Glycolysis = c("ENO1", "ENO2", "ENO3", /* ... */),
  MHC_ClassII = c("CD74", "HLA-DRA", /* ... */),
  Ribosomal_Proteins = c("MRPL1", "MRPL2", /* ... */),
  Coagulation = c("F10", "F11", "F12", /* ... */),
  MHC_ClassI = c("HLA-A", "HLA-B", /* ... */),
  TypeI_IFN_Response = c("ADAR", "CACTIN", /* ... */),
  IFN_Gamma_Response = c("ACTG1", "ACTR2", /* ... */),
  IFN_Beta_Response = c("AIM2", "BST2", /* ... */)
)

# Define cell type groups
CELL_TYPES <- list(
c('CD14+ Mono', 'FCGR3A+ Mono', 'Platelet', 'B', 'CD4+ T', 
          'CD8+ T', 'NK', 'Erythroid precursor', 'DC', 'Neutrophil'))


#' Calculate module scores
#' @param seurat_obj Seurat object
#' @param modules List of gene modules
#' @return Seurat object with module scores
calculate_module_scores <- function(seurat_obj, modules) {
  for (module_name in names(modules)) {
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(modules[[module_name]]),
      name = module_name,
      search = TRUE
    )
  }
  return(seurat_obj)
}

#' Create boxplot for cell type and module
#' @param data Module score data
#' @param cell_type Cell type to plot
#' @param module_name Name of module
#' @return ggplot object
create_module_boxplot <- function(data, cell_type, module_name) {
  # Filter data
  plot_data <- data %>%
    filter(CellType == cell_type) %>%
    select(time, Module) %>%
    mutate(Module = as.numeric(Module))
  
  # Calculate statistics
  stats <- plot_data %>%
    wilcox_test(Module ~ time) %>%
    adjust_pvalue() %>%
    add_significance("p.adj") %>%
    add_xy_position(fun = "max", x = "time")
  
  # Create plot
  ggboxplot(
    plot_data,
    x = "time",
    y = "Module",
    fill = "time",
    xlab = cell_type,
    palette = c("#ef8a62", "#999999")
  ) +
    stat_pvalue_manual(
      stats,
      label = "p.adj.signif",
      tip.length = 0.01
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks.x = element_blank()
    )
}

#' Process module scores for one condition
#' @param seurat_obj Seurat object
#' @param module_name Module name
#' @param condition Condition to analyze
#' @param cell_types Cell types to include
#' @param output_dir Output directory
analyze_condition <- function(seurat_obj, module_name, condition, 
                            cell_types, output_dir) {
  # Extract data
  module_data <- seurat_obj[[c(
    paste0(module_name, "1"),
    "condition",
    "CellType",
    "time"
  )]] %>%
    filter(condition == !!condition)
  
  # Create plots
  plots <- map(cell_types, ~create_module_boxplot(
    module_data,
    .x,
    module_name
  ))
  
  # Save plot
  n_cols <- if(length(cell_types) <= 3) length(cell_types) else 10
  width <- 4 * n_cols
  
  png(
    file.path(output_dir, sprintf("%s_%s_module_scores.png", condition, module_name)),
    width = width * 300,
    height = 4 * 300,
    res = 300,
    type = 'cairo'
  )
  
  print(plot_grid(plotlist = plots, ncol = n_cols))
  
  dev.off()
  
  return(plots)
}

#' Main execution function
#' @param input_path Path to Seurat object
#' @param output_dir Output directory
main <- function(input_path, output_dir) {
  # Load data
  seurat_obj <- readRDS(input_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Add metadata
  seurat_obj <- add_condition_metadata(seurat_obj)
  
  # Calculate module scores
  seurat_obj <- calculate_module_scores(seurat_obj, GENE_MODULES)
  
  # Analyze each condition
  conditions <- c("NS LS", "NS ES", "S")
  
  for (condition in conditions) {
    # Analyze all cell types
    analyze_condition(
      seurat_obj,
      "OXPHOS",
      condition,
      CELL_TYPES$all,
      output_dir
    )
    
    # Analyze lymphoid cells
    analyze_condition(
      seurat_obj,
      "Exhaustion.Score",
      condition,
      CELL_TYPES$lymphoid,
      output_dir
    )
    
    # Analyze cytotoxic cells
    analyze_condition(
      seurat_obj,
      "Cytotoxicity.Score",
      condition,
      CELL_TYPES$cytotoxic,
      output_dir
    )
  }
}

# Example usage:
# main(
#   input_path = "results/seurat_object.rds",
#   output_dir = "results/module_scores"
# )
