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
  library(SeuratDisk)})

#==============================================================================
# Environment variables and helper functions
#==============================================================================

setwd("/scratch/hnatri/PIPAC_spatial/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

#==============================================================================
# Import data and convert
#==============================================================================

cell_seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat.rds")
denoist_seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat_denoIST.rds")

# Seurat v5 assay causes an error
DefaultAssay(cell_seurat_data) <- "RNA"
cell_seurat_data[["RNA"]] <- as(object = cell_seurat_data[["RNA"]], Class = "Assay")
cell_seurat_data <- FindVariableFeatures(cell_seurat_data)
cell_seurat_data[["nucleus_RNA"]] <- as(object = cell_seurat_data[["nucleus_RNA"]], Class = "Assay")
cell_seurat_data <- FindVariableFeatures(cell_seurat_data)

cell_seurat_data[["denoist_RNA"]] <- CreateAssayObject(counts = denoist_seurat_data@assays$denoist_RNA$counts)

cell_seurat_data[["denoist_RNA"]] <- as(object = cell_seurat_data[["denoist_RNA"]], Class = "Assay")
cell_seurat_data <- FindVariableFeatures(cell_seurat_data)

DefaultAssay(cell_seurat_data) <- "denoist_RNA"

cell_seurat_data[["RNA"]] <- NULL
cell_seurat_data[["nucleus_RNA"]] <- NULL
cell_seurat_data[["RNA"]] <- cell_seurat_data[["denoist_RNA"]]
DefaultAssay(cell_seurat_data) <- "RNA"
cell_seurat_data[["denoist_RNA"]] <- NULL

#cell_seurat_data[["RNA"]] <- as(object = cell_seurat_data[["RNA"]], Class = "Assay5")

head(cell_seurat_data@meta.data)
cell_seurat_data@meta.data <- cell_seurat_data@meta.data[,c("cell_id", "Sample")]

cell_seurat_data@assays$RNA@counts
cell_seurat_data@assays$RNA@data

#cell_seurat_data@assays$RNA$data <- NULL
cell_seurat_data@assays$RNA$scale.data <- NULL

cell_seurat_data <- NormalizeData(cell_seurat_data)

#cell_seurat_data[["RNA_filtered"]] <- CreateAssayObject(counts = rna_data)

SaveH5Seurat(cell_seurat_data, filename = "denoist_only_seurat_data.h5Seurat", overwrite = TRUE)
Convert("denoist_only_seurat_data.h5Seurat", dest = "h5ad", overwrite = TRUE)

###

cell_seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat.rds")

# Seurat v5 assay causes an error
DefaultAssay(cell_seurat_data) <- "RNA"
cell_seurat_data[["RNA"]] <- as(object = cell_seurat_data[["RNA"]], Class = "Assay")
cell_seurat_data <- FindVariableFeatures(cell_seurat_data)
cell_seurat_data[["nucleus_RNA"]] <- as(object = cell_seurat_data[["nucleus_RNA"]], Class = "Assay")
cell_seurat_data <- FindVariableFeatures(cell_seurat_data)

DefaultAssay(cell_seurat_data) <- "RNA"

cell_seurat_data[["nucleus_RNA"]] <- NULL

#cell_seurat_data[["RNA"]] <- as(object = cell_seurat_data[["RNA"]], Class = "Assay5")

head(cell_seurat_data@meta.data)
cell_seurat_data@meta.data <- cell_seurat_data@meta.data[,c("cell_id", "Sample")]

cell_seurat_data@assays$RNA@counts
cell_seurat_data@assays$RNA@data

#cell_seurat_data@assays$RNA$data <- NULL
cell_seurat_data@assays$RNA$scale.data <- NULL

#cell_seurat_data[["RNA_filtered"]] <- CreateAssayObject(counts = rna_data)

SaveH5Seurat(cell_seurat_data, filename = "cell_only_seurat_data_test.h5Seurat", overwrite = TRUE)
Convert("cell_only_seurat_data_test.h5Seurat", dest = "h5ad", overwrite = TRUE)


