# canonical_annotation.R
# Annotate cell clusters using canonical marker genes
# Data from Qiu et al., 2021

library(dplyr)

#' Define canonical markers for each cell type
#' @return Named list of marker genes
get_canonical_markers <- function() {
  list(
    "B" = "MS4A1",
    "CD14_Mono" = c("CD14", "LYZ"),
    "CD4_T" = c("IL7R", "CCR7", "CD27"),
    "CD8_T" = "CD8A",
    "DC" = c("FCER1A", "CST3", "CD123", "GZMB"),
    "Erythroid" = c("GYPB", "AHSP"),
    "FCGR3A_Mono" = c("FCGR3A", "MS4A7"),
    "Neutrophil" = c("JAML", "SERPINB"),
    "NK" = c("GNLY", "NKG7"),
    "Platelet" = "PPBP"
  )
}

#' Process marker gene results and annotate clusters
#' @param marker_file Path to marker gene results
#' @param output_dir Output directory
#' @return Annotated cluster information
annotate_clusters <- function(marker_file, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read marker gene results
  markers_df <- read.delim(
    marker_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t"
  )
  
  # Summarize genes by cluster
  cluster_genes <- markers_df %>%
    select(cluster, gene) %>%
    group_by(cluster) %>%
    summarise(genes = paste(gene, collapse = ","))
  
  # Get canonical markers
  canonical_markers <- get_canonical_markers()
  
  # Function to check marker presence
  check_markers <- function(gene_list, markers) {
    all(sapply(markers, function(x) grepl(x, gene_list, fixed = TRUE)))
  }
  
  # Annotate clusters
  cluster_genes <- cluster_genes %>%
    mutate(
      cell_type = case_when(
        check_markers(genes, canonical_markers$B) ~ "B_cell",
        check_markers(genes, canonical_markers$CD14_Mono) ~ "CD14_Monocyte",
        check_markers(genes, canonical_markers$CD4_T) ~ "CD4_T_cell",
        check_markers(genes, canonical_markers$CD8_T) ~ "CD8_T_cell",
        check_markers(genes, canonical_markers$DC) ~ "Dendritic_cell",
        check_markers(genes, canonical_markers$Erythroid) ~ "Erythroid_precursor",
        check_markers(genes, canonical_markers$FCGR3A_Mono) ~ "FCGR3A_Monocyte",
        check_markers(genes, canonical_markers$Neutrophil) ~ "Neutrophil",
        check_markers(genes, canonical_markers$NK) ~ "NK_cell",
        check_markers(genes, canonical_markers$Platelet) ~ "Platelet",
        TRUE ~ "Unknown"
      )
    )
  
  # Save results
  write.table(
    cluster_genes,
    file = file.path(output_dir, "cluster_annotations.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  return(cluster_genes)
}

#' Generate summary of cluster annotations
#' @param annotations Cluster annotation dataframe
#' @param output_dir Output directory
summarize_annotations <- function(annotations, output_dir) {
  # Create summary
  summary <- annotations %>%
    group_by(cell_type) %>%
    summarise(
      n_clusters = n(),
      clusters = paste(cluster, collapse = ",")
    )
  
  # Save summary
  write.table(
    summary,
    file = file.path(output_dir, "annotation_summary.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  return(summary)
}

# Main execution
main <- function(marker_file, output_dir) {
  # Annotate clusters
  annotations <- annotate_clusters(marker_file, output_dir)
  
  # Generate summary
  summary <- summarize_annotations(annotations, output_dir)
  
  # Print results
  print("Annotation complete. Results saved in:")
  print(output_dir)
  
  return(list(
    annotations = annotations,
    summary = summary
  ))
}

# Example usage:
# results <- main(
#   marker_file = "results/markers/marker_genes_mast.txt",
#   output_dir = "results/annotations"
# )
