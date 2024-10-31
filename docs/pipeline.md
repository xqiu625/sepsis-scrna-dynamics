# Single-cell RNA-seq Analysis Pipeline for Sepsis Study

This document outlines the computational pipeline used in the paper "Dynamic changes in human single-cell transcriptional signatures during fatal sepsis" (Qiu et al., 2021).

## Overview

This pipeline analyzes single-cell RNA sequencing data from peripheral blood mononuclear cells (PBMCs) of sepsis patients to understand the molecular dynamics during the critical early hours of sepsis progression. The analysis includes 4 sepsis samples and 2 healthy controls, resulting in 57,133 cells after quality control.

## Data Processing Steps

### 1. Raw Data Processing (Cell Ranger v3.1.0)
```bash
# Cell Ranger processing for sample demultiplexing and counting
cellranger count \
    --id=run_count_sample \
    --fastqs=/path/to/fastqs \
    --sample=sample_id \
    --transcriptome=/path/to/reference
```
Key steps performed by Cell Ranger:
- Sample demultiplexing
- Barcode processing
- 5' UMI counting
- FASTQ alignment using STAR aligner
- Automatic cell barcode determination based on UMI distribution

### 2. Quality Control and Cell Filtering
```R
# Load libraries
library(Seurat)
library(dplyr)

# Read in data
seurat_obj <- Read10X("path/to/cellranger/output")
seurat_obj <- CreateSeuratObject(seurat_obj)

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Apply QC filters
seurat_obj <- subset(seurat_obj,
    subset = nFeature_RNA >= 200 &
            nFeature_RNA <= 6000 &
            nCount_RNA > 1000 &
            percent.mt < 0.2
)
```
Quality control criteria:
- Gene count per cell: 200-6000
- UMI count per cell: >1000
- Mitochondrial gene percentage: <20%

### 3. Data Integration and Batch Correction (Seurat v3)
```R
# Split object by sample
obj.list <- SplitObject(seurat_obj, split.by = "sample")

# Normalize and identify variable features for each dataset
obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x)
})

# Select features for integration
features <- SelectIntegrationFeatures(object.list = obj.list)

# Find integration anchors and integrate data
anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                 anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

# Scale data and run PCA
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined)
combined <- RunPCA(combined)
```

### 4. Downstream Analysis
```R
# Run UMAP
combined <- RunUMAP(combined, dims = 1:30)

# Clustering
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 1.0)

# Save processed object
saveRDS(combined, "sepsis_singlecell_transcriptome.rds")
```

## File Structure
```
sepsis-scrna-dynamics/
├── data/
│   ├── raw/           # Raw FASTQ files
│   ├── cellranger/    # Cell Ranger output
│   └── processed/     # Processed Seurat object
├── src/
│   ├── preprocessing/
│   │   ├── run_cellranger.sh
│   │   └── create_seurat_object.R
│   └── analysis/
└── results/
```

## Dependencies
- Cell Ranger v3.1.0
- R v4.0.0
- Seurat v3.1.5
- See environment.yml for complete list

## Runtime Requirements
- Memory: 100GB recommended
- CPU: 10 cores recommended
- Storage: 100GB minimum

## Dataset Summary
- Sample Groups:
  * Nonsepsis control (n = 2)
  * Sepsis nonsurvivor (n = 2)
  * Sepsis survivor (n = 3)
- Total number of cells after QC: 57,133
- Time points for sepsis samples:
  * T0 (diagnosis)
  * T6 (6 hours post-diagnosis)
- Condition labels:
  * HC: Nonsepsis control
  * NS: Sepsis nonsurvivor
  * S: Sepsis survivor
- Sample labeling:
  * HC1, HC2: Nonsepsis controls
  * NS1_T0, NS1_T6, NS2_T0, NS2_T6: Sepsis nonsurvivors at 0h and 6h
  * S1_T0, S1_T6, S2_T0, S2_T6, S3_T0, S3_T6: Sepsis survivors at 0h and 6h

## Citation
```bibtex
@article{qiu2021dynamic,
  title={Dynamic changes in human single-cell transcriptional signatures during fatal sepsis},
  author={Qiu, Xinru and Li, Jiang and Bonenfant, Jeff and Jaroszewski, Lukasz and Mittal, Aarti and Klein, Walter and Godzik, Adam and Nair, Meera G},
  journal={Journal of Leukocyte Biology},
  volume={110},
  number={6},
  pages={1253--1268},
  year={2021}
}
```
