# singler_annotation.R
# Automated cell type annotation using SingleR with Human Primary Cell Atlas reference
# Data from Qiu et al., 2021

library(SingleCellExperiment)
library(Seurat)
library(celldex)
library(SingleR)
library(dplyr)

#' Run SingleR annotation pipeline
#' @param seurat_obj Seurat object
#' @param output_dir Directory for output files
#' @return List containing annotated metadata and cluster annotations
run_singler_annotation <- function(seurat_obj, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load reference dataset
  ref.data <- HumanPrimaryCellAtlasData(ensembl = FALSE)
  
  # Convert Seurat to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Run SingleR
  pred.hpca <- SingleR(
    test = sce,
    ref = ref.data,
    assay.type.test = 1,
    labels = ref.data$label.main
  )
  
  # Process results
  seurat <- as.Seurat(sce)
  SingleR_metadata <- seurat@meta.data
  
  # Create annotation dataframe
  df <- as.data.frame(pred.hpca) %>%
    dplyr::select(pruned.labels)
  
  # Merge annotations with metadata
  SingleR_metadata <- merge(SingleR_metadata, df, by = 0)
  
  # Save complete metadata
  write.table(
    SingleR_metadata,
    file = file.path(output_dir, "SingleR_metadata.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  # Process cluster-level annotations
  cluster_annotations <- process_cluster_annotations(SingleR_metadata)
  
  # Save cluster annotations
  write.table(
    cluster_annotations,
    file = file.path(output_dir, "singleR_cellType.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(list(
    metadata = SingleR_metadata,
    cluster_annotations = cluster_annotations
  ))
}

#' Process cluster-level annotations
#' @param metadata SingleR metadata
#' @return Dataframe of cluster annotations
process_cluster_annotations <- function(metadata) {
  # Select relevant columns
  cluster_type <- metadata %>% 
    dplyr::select(seurat_clusters, SingleR.labels)
  
  # Get number of clusters
  n_clusters <- n_distinct(cluster_type$seurat_clusters)
  
  # Initialize list to store results
  cluster_results <- list()
  
  # Process each cluster
  for (i in 0:(n_clusters-1)) {
    # Get cell type frequencies for cluster
    cluster_counts <- cluster_type %>%
      dplyr::filter(seurat_clusters == i) %>%
      pull(SingleR.labels) %>%
      table() %>%
      as.data.frame()
    
    # Sort by frequency and get top 2 cell types
    top_types <- cluster_counts %>%
      arrange(desc(Freq)) %>%
      head(2) %>%
      mutate(cluster = i)
    
    # Store results
    cluster_results[[i+1]] <- top_types
  }
  
  # Combine results
  cluster_annotations <- bind_rows(cluster_results) %>%
    rename(
      cellType = Var1,
      cellCount = Freq,
      cluster = cluster
    )
  
  return(cluster_annotations)
}

# Main execution
main <- function() {
  # Set paths
  input_path <- "results/integrated/sepsis_integrated.rds"
  output_dir <- "results/singler"
  
  # Load data
  seurat_obj <- readRDS(input_path)
  
  # Run SingleR annotation
  results <- run_singler_annotation(seurat_obj, output_dir)
  
  # Print summary
  print("SingleR annotation complete. Results saved in:")
  print(output_dir)
}

# Example usage:
# main()
