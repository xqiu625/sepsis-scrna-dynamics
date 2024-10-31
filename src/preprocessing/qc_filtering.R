# qc_filtering.R
# QC filtering for sepsis single-cell RNA-seq data
# Parameters from Qiu et al., 2021

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

#' Read and process 10X data
#' @param data_dir Directory containing Cell Ranger output
#' @param sample_name Name of the sample
#' @return Seurat object with QC metrics
create_seurat_object <- function(data_dir, sample_name) {
  # Read data
  data <- Read10X(data.dir = data_dir)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  return(seurat_obj)
}

#' Filter cells based on QC metrics
#' @param seurat_obj Seurat object
#' @return Filtered Seurat object
filter_cells <- function(seurat_obj) {
  # Apply QC thresholds as specified in the paper
  seurat_filtered <- subset(seurat_obj,
    subset = nFeature_RNA >= 200 &
      nFeature_RNA <= 6000 &
      nCount_RNA > 1000 &
      percent.mt < 0.2
  )
  
  return(seurat_filtered)
}

# Define sample information
sample_info <- list(
  input_dirs = c(
    "cellranger/HC1",
    "cellranger/HC2",
    "cellranger/NS_LS_T0",
    "cellranger/NS_LS_T6",
    "cellranger/NS_ES_T0",
    "cellranger/NS_ES_T6",
    "cellranger/S1_T0",
    "cellranger/S1_T6",
    "cellranger/S2_T0",
    "cellranger/S2_T6",
    "cellranger/S3_T0",
    "cellranger/S3_T6"
  ),
  sample_names = c(
    "HC1",
    "HC2",
    "NS_LS_T0",
    "NS_LS_T6",
    "NS_ES_T0",
    "NS_ES_T6",
    "S1_T0",
    "S1_T6",
    "S2_T0",
    "S2_T6",
    "S3_T0",
    "S3_T6"
  )
)


process_samples <- function(input_dirs, sample_names, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize list for storing QC results
  qc_results <- list()
  filtered_objects <- list()
  
  # Process each sample
  for (i in seq_along(sample_names)) {
    message(paste0("Processing sample: ", sample_names[i]))
    
    # Create and filter Seurat object
    seurat_obj <- create_seurat_object(input_dirs[i], sample_names[i])
    filtered_obj <- filter_cells(seurat_obj)
    
    # Store QC metrics
    qc_results[[i]] <- data.frame(
      Sample = sample_names[i],
      Total_Cells = ncol(seurat_obj),
      Cells_After_QC = ncol(filtered_obj),
      Median_Genes = median(filtered_obj$nFeature_RNA),
      Median_UMI = median(filtered_obj$nCount_RNA),
      Median_Mt_Percent = median(filtered_obj$percent.mt)
    )
    
    # Store filtered object
    filtered_objects[[sample_names[i]]] <- filtered_obj
  }
  
  # Combine QC results
  qc_summary <- do.call(rbind, qc_results)
  write.csv(qc_summary, file = file.path(output_dir, "qc_summary.csv"), row.names = FALSE)
  
  # Return filtered objects
  return(filtered_objects)
}

# Run QC pipeline
filtered_objects <- process_samples(
  input_dirs = sample_info$input_dirs,
  sample_names = sample_info$sample_names,
  output_dir = "results/qc"
)

# Add metadata for condition and time point
add_metadata <- function(seurat_obj, sample_name) {
  seurat_obj$condition <- case_when(
    grepl("^HC", sample_name) ~ "HC",
    grepl("^NS_LS", sample_name) ~ "NS_LS",
    grepl("^NS_ES", sample_name) ~ "NS_ES",
    grepl("^S", sample_name) ~ "S"
  )
  
  seurat_obj$Time <- case_when(
    grepl("^HC", sample_name) ~ "NA",
    grepl("T0$", sample_name) ~ "T0",
    grepl("T6$", sample_name) ~ "T6"
  )
  
  return(seurat_obj)
}

# Add metadata to all objects
filtered_objects <- lapply(names(filtered_objects), function(sample_name) {
  add_metadata(filtered_objects[[sample_name]], sample_name)
})

# Merge all objects
merged_seurat <- merge(
  filtered_objects[[1]], 
  y = filtered_objects[2:length(filtered_objects)],
  add.cell.ids = sample_info$sample_names
)

# Save merged object
saveRDS(merged_seurat, "results/qc/sepsis_merged_filtered.rds")

# Print summary
print(table(merged_seurat$condition, merged_seurat$Time))
