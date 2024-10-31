# cluster_composition.R
# Calculate cell counts and percentages per cluster across samples
# Data from Qiu et al., 2021

library(Seurat)
library(tidyr)
library(dplyr)

#' Calculate cell counts per cluster for each sample
#' @param seurat_obj Seurat object
#' @return Dataframe with counts
get_cluster_counts <- function(seurat_obj) {
  FetchData(seurat_obj, vars = c("ident", "sample")) %>%
    dplyr::count(ident, sample) %>%
    tidyr::spread(ident, n) %>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
}

#' Calculate percentages from counts
#' @param count_data Dataframe with counts
#' @return Dataframe with percentages
calculate_percentages <- function(count_data) {
  # Remove sample column for calculations
  count_matrix <- count_data %>% select(-sample)
  
  # Calculate total cells per sample
  total_cells <- rowSums(count_matrix)
  
  # Calculate percentages
  percentages <- count_matrix %>%
    mutate(across(everything(), ~. / total_cells))
  
  # Add back sample column
  percentages <- bind_cols(
    sample = count_data$sample,
    percentages
  )
  
  return(percentages)
}

#' Combine counts and percentages into single output
#' @param counts Count data
#' @param percentages Percentage data
#' @return Combined dataframe
create_combined_output <- function(counts, percentages) {
  # Convert to matrices for easier manipulation
  count_matrix <- counts %>%
    select(-sample) %>%
    t() %>%
    as.data.frame()
  
  percent_matrix <- percentages %>%
    select(-sample) %>%
    t() %>%
    as.data.frame()
  
  # Add sample names from first dataset
  colnames(count_matrix) <- counts$sample
  colnames(percent_matrix) <- percentages$sample
  
  # Combine counts and percentages
  combined_data <- bind_rows(
    Counts = count_matrix,
    Percentages = percent_matrix,
    .id = "Metric"
  )
  
  return(combined_data)
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
  
  # Calculate counts
  message("Calculating cluster counts...")
  counts <- get_cluster_counts(seurat_obj)
  
  # Calculate percentages
  message("Calculating percentages...")
  percentages <- calculate_percentages(counts)
  
  # Combine results
  message("Combining results...")
  combined_data <- create_combined_output(counts, percentages)
  
  # Save results
  write.table(
    combined_data,
    file = file.path(output_dir, "cluster_composition.txt"),
    sep = '\t',
    quote = FALSE,
    row.names = TRUE
  )
  
  # Save separate files for counts and percentages
  write.table(
    counts,
    file = file.path(output_dir, "cluster_counts.txt"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    percentages,
    file = file.path(output_dir, "cluster_percentages.txt"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Analysis complete. Results saved in: ", output_dir)
  
  return(list(
    counts = counts,
    percentages = percentages,
    combined = combined_data
  ))
}

# Example usage:
# results <- main(
#   input_path = "results/integrated/sepsis_integrated.rds",
#   output_dir = "results/cluster_composition"
# )
