# install_packages.R

# Function to check and install CRAN packages
install_cran_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages)
  }
}

# Function to check and install Bioconductor packages
install_bioc_packages <- function(packages) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    BiocManager::install(new_packages)
  }
}

# CRAN packages
cran_packages <- c(
  "Seurat",
  "dplyr",
  "tidyr",
  "ggplot2",
  "cowplot",
  "pheatmap",
  "RColorBrewer",
  "stringr",
  "rstatix",
  "ggpubr",
  "patchwork"
)

# Bioconductor packages
bioc_packages <- c(
  "SingleCellExperiment",
  "SingleR",
  "MAST",
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Hs.eg.db"
)

# Install packages
message("Installing CRAN packages...")
install_cran_packages(cran_packages)

message("Installing Bioconductor packages...")
install_bioc_packages(bioc_packages)

message("Package installation complete.")
