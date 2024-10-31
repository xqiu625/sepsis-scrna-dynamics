# cd52_cd4t_correlation.R
# Analyze correlation between CD52 and CD4+ T cell activation
# Data from Qiu et al., 2021

library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)

# Define T cell activation genes
T_CELL_GENES <- c(
  "ABL1", "ABL2", "ADA", "ADAM8", "ADAM17", "ADK", "ADRM1", "AGER", 
  "AIF1", "ANKLE1", "ANXA1", "AP3B1", "AP3D1", "APBB1IP", "APC", 
  "ARMC5", "ATP7A", "AZI2", "B2M", "BAD", "BATF", "BAX", "BCL2"
  # ... add more T cell genes as needed
)

# Define condition colors
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

#' Create CD52 correlation plot
#' @param data Processed expression data
#' @return ggplot object
create_correlation_plot <- function(data) {
  ggscatter(
    data,
    x = "CD52",
    y = "T_cell_activation",
    color = "Condition",
    size = 0.1,
    add = "reg.line",
    conf.int = TRUE,
    palette = CONDITION_COLORS
  ) +
    stat_cor(
      aes(color = Condition),
      label.x = 1.5
    ) +
    labs(
      title = "CD4+ T Cell",
      x = "CD52 Expression",
      y = "T Cell Activation Score"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}

#' Process and analyze CD4+ T cell data
#' @param seurat_obj Seurat object
#' @return Processed data frame
process_cd4t_data <- function(seurat_obj) {
  # Subset CD4+ T cells
  cd4t_obj <- subset(seurat_obj, CellType == "CD4+ T")
  
  # Calculate T cell activation score
  cd4t_obj <- AddModuleScore(
    object = cd4t_obj,
    features = list(T_CELL_GENES),
    name = "T_cell_activation",
    search = TRUE
  )
  
  # Extract data
  FetchData(
    cd4t_obj,
    vars = c("T_cell_activation1", "condition", "CD52")
  ) %>%
    as.data.frame() %>%
    filter(CD52 != 0) %>%
    rename(
      T_cell_activation = T_cell_activation1,
      Condition = condition
    ) %>%
    mutate(
      across(c(T_cell_activation, CD52), as.numeric)
    )
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
  message("Processing metadata...")
  seurat_obj <- add_condition_metadata(seurat_obj)
  
  # Process data
  message("Analyzing CD4+ T cells...")
  cd4t_data <- process_cd4t_data(seurat_obj)
  
  # Create plot
  message("Generating correlation plot...")
  plot <- create_correlation_plot(cd4t_data)
  
  # Save plot
  message("Saving results...")
  ggsave(
    file.path(output_dir, "cd52_cd4t_correlation.png"),
    plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # Save processed data
  write_csv(
    cd4t_data,
    file.path(output_dir, "cd52_cd4t_data.csv")
  )
  
  message("Analysis complete")
  
  return(list(
    data = cd4t_data,
    plot = plot
  ))
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/cd52_correlation"
# )
