# violin_plot.R
# Create violin plot for gene expression across conditions
# Data from Qiu et al., 2021

library(Seurat)
library(ggplot2)

#' Create custom violin plot
#' @param seurat_obj Seurat object
#' @param gene Gene to plot
#' @param group_var Variable to group by
#' @return ggplot object
create_violin_plot <- function(seurat_obj, 
                             gene = "HIF1A", 
                             group_var = "condition") {
  # Create base violin plot
  plot <- VlnPlot(
    seurat_obj,
    features = gene,
    sort = 'increasing',
    group.by = group_var,
    adjust = 1,
    pt.size = 0,
    combine = TRUE
  )
  
  # Customize theme
  plot + theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none"
  )
}

#' Save plot with consistent parameters
#' @param plot ggplot object
#' @param output_file Output file path
save_plot <- function(plot, output_file) {
  png(
    output_file,
    width = 4 * 300,
    height = 4 * 300,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  print(plot)
  dev.off()
}

#' Main execution function
#' @param seurat_obj Seurat object
#' @param output_dir Output directory
main <- function(seurat_obj, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create plot
  violin_plot <- create_violin_plot(seurat_obj)
  
  # Save plot
  save_plot(
    violin_plot,
    file.path(output_dir, "hif1a_violin.png")
  )
  
  return(violin_plot)
}

# Example usage:
# plot <- main(
#   seurat_obj = seurat_subset,
#   output_dir = "results/figures"
# )


