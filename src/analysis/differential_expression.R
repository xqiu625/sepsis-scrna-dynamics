# differential_expression.R

library(Seurat)
library(tidyverse)
library(MAST)
library(scCATCH)


find_marker_genes <- function(seurat_obj, output_dir) {
  # Set RNA as default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find markers using MAST
  markers <- FindAllMarkers(
    object = seurat_obj,
    assay = 'RNA',
    only.pos = TRUE,
    test.use = 'MAST'
  )
  
  # Save marker results
  write.table(
    markers,
    file = file.path(output_dir, "marker_genes_mast.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = '\t'
  )
  
  return(list(
    seurat_obj = seurat_obj,
    markers = markers
  ))
}
