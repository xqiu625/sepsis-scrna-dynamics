# module_dotplot_analysis.R
# Calculate module scores and create dot plots
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)

# Define gene modules
GENE_MODULES <- list(
  OXPHOS = c(
    "NDUFA10", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
    # ... [rest of OXPHOS genes]
  ),
  
  Glycolysis = c(
    "ENO1", "ENO2", "ENO3", "ENTPD5", "GAPDHS",
    # ... [rest of Glycolysis genes]
  )
)

# Define sample type order
SAMPLE_ORDER <- c(
  "S T6", "S T0",
  "NS ES T6", "NS ES T0",
  "NS LS T6", "NS LS T0",
  "HC"
)

  
  # Set factor levels for proper ordering
  metadata$sample_type <- factor(
    metadata$sample_type,
    levels = SAMPLE_ORDER
  )
  
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}

#' Calculate module scores
#' @param seurat_obj Seurat object
#' @param gene_modules List of gene modules
#' @return Seurat object with module scores
calculate_module_scores <- function(seurat_obj, gene_modules) {
  for (module_name in names(gene_modules)) {
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(gene_modules[[module_name]]),
      name = module_name,
      search = TRUE
    )
    
    # Rename score column
    old_name <- paste0(module_name, "1")
    colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old_name] <- module_name
  }
  
  return(seurat_obj)
}

#' Create module score dot plot
#' @param seurat_obj Seurat object
#' @param cell_type Cell type to analyze
#' @param output_file Output file path
create_dotplot <- function(seurat_obj, cell_type, output_file) {
  # Subset data
  subset_obj <- subset(seurat_obj, idents = cell_type)
  
  # Create plot
  plot <- DotPlot(
    object = subset_obj,
    features = names(GENE_MODULES),
    cols = c("lightgrey", "blue"),
    col.min = -1.0,
    col.max = 1.0,
    dot.min = 0,
    dot.scale = 6,
    group.by = "sample_type"
  ) +
    RotatedAxis() +
    theme(
      axis.title = element_blank()
    )
  
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
  seurat_obj <- add_metadata(seurat_obj)
  seurat_obj$CellType <- Idents(seurat_obj)
  
  # Calculate module scores
  message("Calculating module scores...")
  seurat_obj <- calculate_module_scores(seurat_obj, GENE_MODULES)
  
  # Create plots for each cell type
  cell_types <- c("CD14+ Mono", "FCGR3A+ Mono")
  
  plots <- map(cell_types, function(cell_type) {
    message("Creating plot for ", cell_type)
    output_file <- file.path(
      output_dir,
      paste0(gsub(" ", "_", cell_type), "_ATP_module.png")
    )
    
    create_dotplot(seurat_obj, cell_type, output_file)
  })
  
  names(plots) <- cell_types
  return(plots)
}

# Example usage:
# results <- main(
#   input_path = "data/sepsis_singlecell_transcriptome.rds",
#   output_dir = "results/module_scores"
# )
