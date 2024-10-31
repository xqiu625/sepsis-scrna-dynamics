# go_term_analysis.R
# Analyze and visualize GO terms using REVIGO and heatmaps
# Data from Qiu et al., 2021

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(pheatmap)
library(RColorBrewer)

# Define color scheme
COLOR_BREAKS <- c(seq(-11, -5.160504e-16, by=1), seq(5.160504e-18, 11, by=1))
COLOR_PALETTE <- c(
  colorRampPalette(c("#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", 
                     "#a6bddb", "#d0d1e6", "#ece7f2", "#fff7fb"))(length(COLOR_BREAKS)/2),
  colorRampPalette(c("#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59", 
                     "#ef6548", "#d7301f", "#b30000", "#7f0000"))(length(COLOR_BREAKS)/2)
)

#' Process GO analysis results using REVIGO
#' @param go_file Path to GO analysis results
#' @param top_n Number of top terms to keep
#' @param negative Whether to make scores negative
#' @return Processed GO terms data frame
process_go_terms <- function(go_file, top_n = 10, negative = FALSE) {
  # Read and process GO results
  go_analysis <- read.delim(go_file)
  
  # Calculate similarity matrix
  sim_matrix <- calculateSimMatrix(
    go_analysis$ID,
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel"
  )
  
  # Calculate scores and reduce terms
  scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
  reduced_terms <- reduceSimMatrix(
    sim_matrix,
    scores,
    threshold = 0.4,
    orgdb = "org.Hs.eg.db"
  )
  
  # Process results
  df <- as.data.frame(reduced_terms) %>%
    select(parentTerm, score) %>%
    mutate(
      parentTerm = paste0("GO: ", parentTerm)
    ) %>%
    arrange(desc(score)) %>%
    distinct(parentTerm, .keep_all = TRUE) %>%
    slice_head(n = top_n)
  
  # Make scores negative if requested
  if (negative) {
    df <- df %>%
      arrange(score) %>%
      mutate(score = -score)
  }
  
  return(df)
}

#' Create GO term heatmap
#' @param data Matrix of GO term scores
#' @param output_file Output file path
#' @param width Plot width
create_heatmap <- function(data, output_file, width = 0.5) {
  # Calculate height based on number of terms
  height <- dim(data)[1] * 0.3
  
  # Create plot
  png(
    output_file,
    width = 300 * width,
    height = 300 * height,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  
  pheatmap(
    data,
    border_color = NA,
    cellwidth = 20,
    cellheight = 20,
    color = COLOR_PALETTE,
    breaks = COLOR_BREAKS,
    legend = FALSE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
  
  dev.off()
}

#' Create scale legend
#' @param output_file Output file path
create_scale_legend <- function(output_file) {
  scale <- matrix(-11:11, nrow = 23, ncol = 1)
  
  png(
    output_file,
    width = 1500,
    height = 300,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  
  pheatmap(
    scale,
    border_color = NA,
    cellwidth = 20,
    cellheight = 20,
    color = COLOR_PALETTE,
    breaks = COLOR_BREAKS,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
  
  dev.off()
}

#' Process comparison between two conditions
#' @param condition1_file Path to first condition's GO results
#' @param condition2_file Path to second condition's GO results
#' @param output_prefix Prefix for output files
#' @param output_dir Output directory
process_comparison <- function(condition1_file, condition2_file, 
                             output_prefix, output_dir) {
  # Process both conditions
  df1 <- process_go_terms(condition1_file, negative = TRUE)
  df2 <- process_go_terms(condition2_file)
  
  # Combine results
  combined_df <- rbind(df2, df1)
  
  # Save term list
  write.csv(
    combined_df,
    file.path(output_dir, paste0(output_prefix, "_terms.csv"))
  )
  
  # Prepare matrix for heatmap
  heatmap_matrix <- combined_df %>%
    column_to_rownames("parentTerm") %>%
    as.matrix()
  
  # Create heatmap
  create_heatmap(
    heatmap_matrix,
    file.path(output_dir, paste0(output_prefix, "_heatmap.png"))
  )
  
  return(list(
    data = combined_df,
    matrix = heatmap_matrix
  ))
}

#' Main execution function
#' @param sepsis_hc_dir Directory with Sepsis vs HC results
#' @param ns_s_dir Directory with NS vs S results
#' @param output_dir Output directory
main <- function(sepsis_hc_dir, ns_s_dir, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process Sepsis vs HC comparison
  sepsis_hc_results <- process_comparison(
    file.path(sepsis_hc_dir, "Erythroblast_HC_bp.txt"),
    file.path(sepsis_hc_dir, "Erythroblast_Sepsis_bp.txt"),
    "sepsis_vs_hc",
    output_dir
  )
  
  # Process NS vs S comparison
  ns_s_results <- process_comparison(
    file.path(ns_s_dir, "Erythroblast_Sepsis_S_bp.txt"),
    file.path(ns_s_dir, "Erythroblast_Sepsis_NS_bp.txt"),
    "ns_vs_s",
    output_dir
  )
  
  # Create scale legend
  create_scale_legend(file.path(output_dir, "scale_legend.png"))
  
  return(list(
    sepsis_hc = sepsis_hc_results,
    ns_s = ns_s_results
  ))
}

# Example usage:
# results <- main(
#   sepsis_hc_dir = "results/go_analysis/sepsis_vs_hc",
#   ns_s_dir = "results/go_analysis/ns_vs_s",
#   output_dir = "results/go_visualization"
# )
