# singler_annotation.R
# Automatic cell type annotation using SingleR
# Data from Qiu et al., 2021

library(SingleR)
library(Seurat)
library(tidyverse)
library(celldex)

# Load reference dataset
ref_data <- HumanPrimaryCellAtlasData()

# Function to perform SingleR annotation
run_singler_annotation <- function(seurat_obj, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract expression matrix from Seurat object
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  
  # Run SingleR
  pred <- SingleR(
    test = expr_matrix,
    ref = ref_data,
    labels = ref_data$label.main,
    method = "cluster",
    clusters = seurat_obj$seurat_clusters
  )
  
  # Add SingleR labels to Seurat object
  seurat_obj$singler_labels <- pred$labels[match(seurat_obj$seurat_clusters,
                                                rownames(pred))]
  
  # Generate label comparison plot
  pdf(file.path(output_dir, "singler_umap.pdf"))
  print(DimPlot(seurat_obj, 
                reduction = "umap", 
                group.by = "singler_labels", 
                label = TRUE))
  dev.off()
  
  # Save annotation results
  write.csv(
    data.frame(
      cluster = rownames(pred),
      singler_label = pred$labels,
      scores = pred$scores
    ),
    file = file.path(output_dir, "singler_annotations.csv"),
    row.names = FALSE
  )
  
  return(seurat_obj)
}

# Main execution
main <- function() {
  # Set paths
  input_path <- "results/integrated/labeled_clustered_integrated_sepsis.rds"
  output_dir <- "results/singler"
  
  # Load integrated Seurat object
  seurat_obj <- readRDS(input_path)
  
  # Run SingleR annotation
  seurat_obj <- run_singler_annotation(seurat_obj, output_dir)
  
  # Save annotated object
  saveRDS(seurat_obj, 
          file = file.path(output_dir, "singler_annotated_sepsis.rds"))
  
  # Compare manual and SingleR annotations
  comparison_df <- data.frame(
    manual = seurat_obj$CellType,
    singler = seurat_obj$singler_labels
  )
  
  write.csv(table(comparison_df), 
            file = file.path(output_dir, "annotation_comparison.csv"))
}

# Run pipeline
main()
