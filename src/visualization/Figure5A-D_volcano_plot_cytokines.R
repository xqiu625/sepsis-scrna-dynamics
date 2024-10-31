# volcano_plot_cytokines.R
# Generate volcano plots highlighting cytokine expression in monocyte subsets
# Data from Qiu et al., 2021

library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(Seurat)

#' Process differential expression results
#' @param deg_data Differential expression results
#' @param control_group Name of control group
#' @return Processed data frame
process_deg_results <- function(deg_data, control_group) {
  deg_data %>%
    mutate(
      across(c(avg_logFC, p_val_adj, p_val), as.numeric),
      # Flip sign if control group
      avg_logFC = ifelse(cluster == control_group, -avg_logFC, avg_logFC),
      # Handle zero p-values
      p_val_adj = ifelse(p_val_adj == 0, 5e-324, p_val_adj)
    ) %>%
    select(gene, p_val_adj, avg_logFC)
}

#' Create enhanced volcano plot
#' @param data Processed DEG data
#' @param cytokine_genes List of cytokine genes to highlight
#' @param output_file Output file path
create_volcano_plot <- function(data, cytokine_genes, output_file) {
  plot <- EnhancedVolcano(
    data,
    lab = data$gene,
    x = 'avg_logFC',
    y = 'p_val_adj',
    selectLab = cytokine_genes,
    xlim = c(min(data$avg_logFC) - 0.5, max(data$avg_logFC) + 0.5),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'adjusted P value'),
    title = "",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.8,
    pointSize = 4.0,
    labSize = 5.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 1,
    caption = "FC cutoff, 0.5; adjusted P value, 0.05",
    legendPosition = 'none',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black'
  )
  
  # Save plot
  png(
    output_file,
    width = 12 * 300,
    height = 6 * 300,
    units = "px",
    res = 300,
    type = 'cairo'
  )
  print(plot)
  dev.off()
  
  return(plot)
}

#' Load and process cytokine database
#' @param file_path Path to cytokine database
#' @return Processed cytokine gene list
load_cytokine_database <- function(file_path) {
  read.csv(file_path) %>%
    rename(hgnc_symbol = 1) %>%
    mutate(
      hgnc_symbol = case_when(
        hgnc_symbol == 'CTF2' ~ 'CEBPZ',
        hgnc_symbol == 'FASL' ~ 'FASLG',
        hgnc_symbol == 'GPI1' ~ 'PIGQ',
        hgnc_symbol == 'KITL' ~ 'KITLG',
        TRUE ~ hgnc_symbol
      )
    )
}

#' Main execution function
#' @param deg_dir Directory containing DEG results
#' @param cytokine_db_path Path to cytokine database
#' @param output_dir Output directory
main <- function(deg_dir, cytokine_db_path, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load cytokine database
  cytokine_db <- load_cytokine_database(cytokine_db_path)
  
  # Analysis configurations
  analyses <- list(
    list(
      name = "CD14_Mono_Sepsis_vs_HC",
      file = "CD14+ Mono_Sepsis_vs_HC.txt",
      control = "HC"
    ),
    list(
      name = "FCGR3A_Mono_Sepsis_vs_HC",
      file = "FCGR3A+ Mono_Sepsis_vs_HC.txt",
      control = "HC"
    ),
    list(
      name = "CD14_Mono_NS_vs_S",
      file = "CD14+ Mono_NS_vs_S.txt",
      control = "Sepsis_S"
    ),
    list(
      name = "FCGR3A_Mono_NS_vs_S",
      file = "FCGR3A+ Mono_NS_vs_S.txt",
      control = "Sepsis_S"
    )
  )
  
  # Process each analysis
  results <- lapply(analyses, function(analysis) {
    # Read and process DEG data
    deg_data <- read.delim2(file.path(deg_dir, analysis$file)) %>%
      process_deg_results(analysis$control)
    
    # Find overlapping cytokine genes
    cytokine_genes <- intersect(deg_data$gene, cytokine_db$hgnc_symbol)
    
    # Create volcano plot
    plot <- create_volcano_plot(
      deg_data,
      cytokine_genes,
      file.path(output_dir, paste0(analysis$name, "_volcano.png"))
    )
    
    return(list(
      data = deg_data,
      cytokines = cytokine_genes,
      plot = plot
    ))
  })
  
  return(results)
}

# Example usage:
# results <- main(
#   deg_dir = "results/differential_expression",
#   cytokine_db_path = "data/cytokine_database.csv",
#   output_dir = "results/volcano_plots"
# )
