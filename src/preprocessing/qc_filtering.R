# qc_filtering.R
# QC filtering for sepsis single-cell RNA-seq data
# Parameters from Qiu et al., 2021

library(Seurat)
library(tidyverse)

#' Create Seurat objects from raw data
#' @param data_dir Base directory containing filtered matrix folders
#' @return List of Seurat objects
create_seurat_objects <- function(data_dir) {
  # Sample names
  samples <- c(
    "HC1", "HC2",
    "NS1_T0", "NS1_T6",
    "NS2_T0", "NS2_T6",
    "S1_T0", "S1_T6",
    "S2_T0", "S2_T6",
    "S3_T0", "S3_T6"
  )
  
  # Create Seurat objects
  seurat_objects <- list()
  for (sample in samples) {
    data_path <- file.path(data_dir, paste0(sample, "_filtered_feature_bc_matrix"))
    seurat_data <- Read10X(data.dir = data_path)
    seurat_obj <- CreateSeuratObject(
      counts = seurat_data,
      min.features = 100,
      project = sample
    )
    seurat_objects[[sample]] <- seurat_obj
  }
  
  return(seurat_objects)
}

#' Merge Seurat objects and add metadata
#' @param seurat_objects List of Seurat objects
#' @return Merged Seurat object with metadata
merge_and_add_metadata <- function(seurat_objects) {
  # Merge objects
  merged_seurat <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[2:length(seurat_objects)],
    add.cell.ids = names(seurat_objects)
  )
  
  # Add metadata
  merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / 
    log10(merged_seurat$nCount_RNA)
  
  # Calculate mitochondrial ratio
  merged_seurat$mitoRatio <- PercentageFeatureSet(
    object = merged_seurat,
    pattern = "^MT-"
  ) / 100
  
  # Add sample information
  merged_seurat$sample <- NA
  for (sample in names(seurat_objects)) {
    merged_seurat$sample[which(str_detect(
      colnames(merged_seurat),
      paste0("^", sample)
    ))] <- sample
  }
  
  return(merged_seurat)
}

#' Generate QC plots
#' @param seurat_obj Seurat object
#' @param output_dir Directory for output plots
plot_qc_metrics <- function(seurat_obj, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata <- seurat_obj@meta.data
  
  # Cell counts per sample
  p1 <- metadata %>% 
    ggplot(aes(x = sample, fill = sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggtitle("Number of Cells per Sample")
  ggsave(file.path(output_dir, "cell_counts.png"), p1, width = 8, height = 6)
  
  # UMI counts
  p2 <- metadata %>% 
    ggplot(aes(x = nCount_RNA, color = sample, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ggtitle("UMI Counts Distribution")
  ggsave(file.path(output_dir, "umi_distribution.png"), p2, width = 8, height = 6)
  
  # Gene counts
  p3 <- metadata %>% 
    ggplot(aes(x = nFeature_RNA, color = sample, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ggtitle("Genes per Cell Distribution")
  ggsave(file.path(output_dir, "genes_per_cell.png"), p3, width = 8, height = 6)
  
  # Mitochondrial ratio
  p4 <- metadata %>% 
    ggplot(aes(x = mitoRatio, color = sample, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    ggtitle("Mitochondrial Content")
  ggsave(file.path(output_dir, "mito_content.png"), p4, width = 8, height = 6)
}

#' Filter cells based on QC metrics
#' @param seurat_obj Seurat object
#' @return Filtered Seurat object
filter_cells <- function(seurat_obj) {
  filtered_obj <- subset(
    x = seurat_obj,
    subset = nFeature_RNA >= 200 &
      nFeature_RNA <= 6000 &
      mitoRatio < 0.20 &
      nCount_RNA > 1000
  )
  
  return(filtered_obj)
}

#' Main execution function
#' @param data_dir Input directory containing filtered matrices
#' @param output_dir Output directory for results
main <- function(data_dir, output_dir) {
  # Create Seurat objects
  seurat_objects <- create_seurat_objects(data_dir)
  
  # Merge and add metadata
  merged_seurat <- merge_and_add_metadata(seurat_objects)
  
  # Generate QC plots before filtering
  plot_qc_metrics(merged_seurat, file.path(output_dir, "pre_qc"))
  
  # Filter cells
  filtered_seurat <- filter_cells(merged_seurat)
  
  # Generate QC plots after filtering
  plot_qc_metrics(filtered_seurat, file.path(output_dir, "post_qc"))
  
  # Save filtered data
  saveRDS(filtered_seurat, file.path(output_dir, "filtered_seurat.rds"))
  
  # Print summary
  print(table(filtered_seurat$sample))
}

# Example usage:
# main(
#   data_dir = "data/filtered_matrices",
#   output_dir = "results/preprocessing"
# )
