# differential_expression_viz.R
# Visualize differential expression between T0 and T6 timepoints
# Data from Qiu et al., 2021

library(tidyverse)
library(ggrepel)
library(RColorBrewer)

# Constants
CELL_TYPES <- c("B", "CD4+ T", "CD8+ T")
GENES_TO_LABEL <- c("CCR2", "CCR6", "CD3D", "CD3E", "CD3G", "CD48", 
                    "CD52", "CD80", "ITGB7", "SELL", "SLAMF6", "THY1")
GROUP_COLORS <- c(
  "0" = "#b3e2cd",
  "1" = "#fdcdac",
  "2" = "#cbd5e8",
  "3" = "#f4cae4",
  "4" = "#e6f5c9"
)

#' Process differential expression results
#' @param file_path Path to DEG results
#' @param cluster_name Cluster name
#' @param position Cluster position
#' @return Processed dataframe
process_deg_file <- function(file_path, cluster_name, position) {
  read.delim2(file_path) %>%
    mutate(
      across(c(avg_logFC, p_val_adj), as.numeric),
      threshold = if_else(cluster == "T6", "Up", "Down"),
      avg_logFC = if_else(cluster == "T6", avg_logFC, -avg_logFC),
      cluster = cluster_name,
      cluster_position = position
    ) %>%
    filter(p_val_adj < 0.05)
}

#' Create background position data
#' @param marker_data Combined marker data
#' @return Background position dataframe
create_background_data <- function(marker_data) {
  marker_data %>%
    group_by(cluster_position) %>%
    summarise(
      Min = min(avg_logFC) - 0.2,
      Max = max(avg_logFC) + 0.2
    ) %>%
    mutate(
      cluster_position = as.numeric(as.vector(cluster_position)),
      start = cluster_position - 0.4,
      end = cluster_position + 0.4
    )
}

#' Create cluster bar position data
#' @param background_data Background position data
#' @return Cluster bar position dataframe
create_cluster_bars <- function(background_data) {
  background_data %>%
    mutate(
      start = cluster_position - 0.5,
      end = cluster_position + 0.5,
      cluster_position = factor(cluster_position, 
                              levels = 0:max(as.vector(cluster_position))),
      group = c("S1", "S2", "S3", "NS LS", "NS ES")
    )
}

#' Create differential expression plot
#' @param marker_data Combined marker data
#' @param background_data Background position data
#' @param cluster_bars Cluster bar data
#' @param labeled_genes Genes to label
#' @param show_legend Whether to show legend
#' @return ggplot object
create_deg_plot <- function(marker_data, background_data, cluster_bars, 
                           labeled_genes, show_legend = FALSE) {
  # Create base plot
  p <- ggplot() +
    # Add background shading
    geom_rect(
      data = background_data,
      aes(xmin = start, xmax = end, ymin = Min, ymax = Max),
      fill = "#525252",
      alpha = 0.1
    ) +
    # Add gene points
    geom_jitter(
      data = marker_data,
      aes(x = cluster_position, y = avg_logFC, color = "orange"),
      size = 1,
      position = position_jitter(seed = 1)
    ) +
    # Set x-axis limits
    scale_x_continuous(
      limits = c(-0.5, max(marker_data$cluster_position) + 0.5),
      breaks = 0:max(marker_data$cluster_position)
    ) +
    # Add gene labels
    geom_label_repel(
      data = filter(marker_data, gene %in% labeled_genes),
      aes(x = cluster_position, y = avg_logFC, label = gene),
      size = 4,
      fontface = "bold",
      segment.curvature = -0.6,
      min.segment.length = 0,
      segment.size = 0.5,
      segment.color = "grey50",
      direction = "y",
      hjust = 1
    ) +
    # Add cluster bars
    geom_rect(
      data = cluster_bars,
      aes(xmin = start, xmax = end, 
          ymin = -0.1, ymax = 0.1,
          fill = cluster_position),
      color = "black",
      alpha = 1,
      show.legend = FALSE
    ) +
    # Add cluster labels
    geom_text(
      data = cluster_bars,
      aes(x = as.numeric(cluster_position) - 1, y = 0, label = group),
      size = 4
    ) +
    scale_fill_manual(values = GROUP_COLORS) +
    theme_bw()
  
  # Add theme elements
  if (!show_legend) {
    p <- p + theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.line.y = element_line(colour = "black")
    )
  }
  
  return(p)
}

#' Process and visualize DEG results for one cell type
#' @param cell_type Cell type to analyze
#' @param input_dirs Directories containing results
#' @param output_dir Output directory
process_cell_type <- function(cell_type, input_dirs, output_dir) {
  # Process survivor results
  s_results <- map2_dfr(
    list(
      file.path(input_dirs$s_dir, paste0(cell_type, "_S1_T6_vs_T0.txt")),
      file.path(input_dirs$s_dir, paste0(cell_type, "_S2_T6_vs_T0.txt")),
      file.path(input_dirs$s_dir, paste0(cell_type, "_S3_T6_vs_T0.txt"))
    ),
    list("S1", "S2", "S3"),
    0:2,
    process_deg_file
  )
  
  # Process non-survivor results
  ns_results <- map2_dfr(
    list(
      file.path(input_dirs$ns_dir, paste0(cell_type, "_NS1_T6_vs_T0.txt")),
      file.path(input_dirs$ns_dir, paste0(cell_type, "_NS2_T6_vs_T0.txt"))
    ),
    list("NS LS", "NS ES"),
    3:4,
    process_deg_file
  )
  
  # Combine results
  combined_results <- bind_rows(s_results, ns_results)
  
  # Create visualization data
  background_data <- create_background_data(combined_results)
  cluster_bars <- create_cluster_bars(background_data)
  
  # Create and save plot
  plot <- create_deg_plot(
    combined_results,
    background_data,
    cluster_bars,
    GENES_TO_LABEL
  )
  
  ggsave(
    file.path(output_dir, paste0(cell_type, "_deg_volcano.png")),
    plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  return(list(
    data = combined_results,
    plot = plot
  ))
}

#' Main execution function
#' @param s_dir Directory with survivor results
#' @param ns_dir Directory with non-survivor results
#' @param output_dir Output directory
main <- function(s_dir, ns_dir, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each cell type
  results <- map(
    CELL_TYPES,
    process_cell_type,
    list(s_dir = s_dir, ns_dir = ns_dir),
    output_dir
  )
  
  names(results) <- CELL_TYPES
  return(results)
}

# Example usage:
# results <- main(
#   s_dir = "results/differential_expression/survivor",
#   ns_dir = "results/differential_expression/non_survivor",
#   output_dir = "results/figures"
# )
