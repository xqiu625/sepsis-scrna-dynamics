# scCATCH.R
# Marker gene identification and cell type annotation
# Data from Qiu et al., 2021

library(tidyverse)
library(scCATCH)

markers <- read.delim2("marker_genes_mast.txt", header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")

annotate_cell_types <- function(markers, output_dir) {
  # Run scCATCH annotation
  cell_types <- scCATCH(
    object = markers,
    species = 'Human',
    cancer = NULL,
    tissue = 'Blood'
  )
  
  # Save annotations
  write.table(
    cell_types,
    file = file.path(output_dir, "cell_type_annotations.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(cell_types)
}

