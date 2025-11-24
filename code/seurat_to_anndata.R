#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 09/04/2025
# Description: Convert Seurat to AnnData with SeuratDisk
#==============================================================================#

#==============================================================================
# Import packages
#==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(anndata)})

#==============================================================================
# Environment variables and helper functions
#==============================================================================

setwd("/home/hnatri/PIPAC_spatial/code/RSC_latest_EDM_2025-08-06/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

#==============================================================================
# Import data and convert
#==============================================================================

cell_seurat_data <- readRDS("/scratch/hnatri/PIPAC/myeloid_seurat_data.rds")
#cell_seurat_immune <- readRDS("/scratch/hnatri/PIPAC/cell_immune_subset.rds")
#cell_seurat_nonimmune <- readRDS("/scratch/hnatri/PIPAC/cell_nonimmune_subset.rds")

#cell_seurat_data <- cell_seurat_immune

# Seurat v5 assay causes an error
DefaultAssay(cell_seurat_data) <- "RNA"
cell_seurat_data[["RNA"]] <- as(object = cell_seurat_data[["RNA"]], Class = "Assay")
cell_seurat_data@assays$RNA$data <- NULL
cell_seurat_data@assays$RNA$scale.data <- NULL

#cell_seurat_data@assays$RNA <- NULL
cell_seurat_data@assays$nucleus_RNA <- NULL
cell_seurat_data@assays$denoist_RNA <- NULL

cell_seurat_data@reductions$umap <- NULL
cell_seurat_data@reductions$pca <- NULL
cell_seurat_data@reductions$sp <- NULL

# NAs in metadata cause a problem
cell_seurat_data@meta.data <- cell_seurat_data@meta.data[,c("cell_id", "Sample")] #, "Annotation"

SaveH5Seurat(cell_seurat_data, filename = "/scratch/hnatri/PIPAC/RSC_latest_EDM_2025-08-06/myeloid.h5Seurat", overwrite = TRUE)
Convert("/scratch/hnatri/PIPAC/RSC_latest_EDM_2025-08-06/myeloid.h5Seurat", dest = "h5ad", overwrite = TRUE)

