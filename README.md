# sepsis-scrna-dynamics

[![DOI](https://zenodo.org/badge/DOI/10.1002/JLB.5MA0721-825R.svg)](https://doi.org/10.1002/JLB.5MA0721-825R)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the analysis code for the paper "Dynamic changes in human single-cell transcriptional signatures during fatal sepsis" (Qiu et al., 2021, Journal of Leukocyte Biology).

## Overview

This project analyzes single-cell RNA sequencing data from peripheral blood mononuclear cells (PBMCs) of sepsis patients to understand the molecular dynamics during the critical early hours of sepsis progression. The analysis includes:

- Processing of 10X Genomics single-cell RNA-seq data
- Cell type identification and annotation
- Differential expression analysis between survivors and non-survivors
- Temporal trajectory analysis of immune cell populations
- Investigation of platelet and erythroid precursor responses
- Analysis of monocyte transcriptional changes
- CD52 expression analysis in lymphocytes

## Data Availability

The raw data have been deposited in the Gene Expression Omnibus (GEO) under accession number [GSE167363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363).

## Prerequisites

### Software Requirements
- R >= 4.0.0
- Cell Ranger >= 3.1.0

### R Packages
```R
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
BiocManager::install(c("SingleR", "MAST", "clusterProfiler"))
```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/sepsis-scrna-dynamics.git
cd sepsis-scrna-dynamics

# Set up R environment
conda env create -f environment.yml
conda activate sepsis-scrna

# Install R dependencies
Rscript src/utils/install_packages.R
```

## Usage

### 1. Data Preprocessing
```bash
# Process raw 10X data
bash src/preprocessing/cellranger_processing.sh

# Perform QC and filtering
Rscript src/preprocessing/qc_filtering.R
```

### 2. Main Analysis
```R
# Run Seurat analysis pipeline
Rscript src/analysis/seurat_integration.R
Rscript src/analysis/cell_type_annotation.R
Rscript src/analysis/differential_expression.R
```

### 3. Generate Figures
```R
# Create visualization plots
Rscript src/visualization/umap_plots.R
Rscript src/visualization/heatmaps.R
```

## Project Structure


For detailed directory structure, see [docs/pipeline.md](docs/pipeline.md).

## Results

The analysis reveals:
- Dynamic changes in immune cell populations within 6 hours of sepsis diagnosis
- Platelet and erythroid precursor responses as drivers of fatal sepsis
- Transcriptional signatures shared with severe COVID-19
- CD52 as a potential prognostic biomarker
- Hypoxic stress as a driving factor in immune dysfunction


## Citation

If you use this code in your research, please cite:

```bibtex
@article{qiu2021dynamic,
  title={Dynamic changes in human single-cell transcriptional signatures during fatal sepsis},
  author={Qiu, Xinru and Li, Jiang and Bonenfant, Jeff and Jaroszewski, Lukasz and Mittal, Aarti and Klein, Walter and Godzik, Adam and Nair, Meera G},
  journal={Journal of Leukocyte Biology},
  volume={110},
  number={6},
  pages={1253--1268},
  year={2021},
  publisher={Wiley Online Library}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Meera Nair** - meera.nair@medsch.ucr.edu
- **Adam Godzik** - adam.godzik@medsch.ucr.edu

## Acknowledgments

This research was supported by:
- UCR School of Medicine
- Dean Innovation Fund
- National Institutes of Health (NIAID, R21AI37830, and R01AI153195)
