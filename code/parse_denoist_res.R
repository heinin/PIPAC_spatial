#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 08/23/2025
# Description: Denoising cellular transcripts
#==============================================================================#

#==============================================================================
# Import packages
#==============================================================================

suppressPackageStartupMessages({
  library(workflowr)
  library(arrow)
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(googlesheets4)
  library(workflowr)})

#==============================================================================
# Environment variables and helper functions
#==============================================================================

setwd("/home/hnatri/PIPAC_spatial/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

source("/home/hnatri/PIPAC_spatial/code/PIPAC_colors_themes.R")
source("/home/hnatri/PIPAC_spatial/code/plot_functions.R")

#==============================================================================
# Import data
#==============================================================================

cell_seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_metadata.rds")

# DenoIST results
res_files <- list.files("/scratch/hnatri/DenoIST/TMA_outputs/")

# Adjusted counts
denoist_counts <- lapply(res_files, function(tma){
  message(tma)
  res <- readRDS(paste0("/scratch/hnatri/DenoIST/TMA_outputs/", tma))
  res <- as.data.frame(res[["adjusted_counts"]])
  
  return(res)
})

lapply(denoist_counts, dim)

#denoist_counts_df <- cbind(denoist_counts)

#tma1 <- readRDS("/scratch/hnatri/DenoIST/TMA_outputs/MR_PIPAC-TMA1")
#tma1$params[[1]]

#==============================================================================
# Add to Seurat
#==============================================================================

seurat_list <- SplitObject(cell_seurat_data, split.by = "TMA")
names(seurat_list) <- unique(cell_seurat_data$TMA)

setdiff(colnames(seurat_list[[1]]), colnames(denoist_counts[[1]]))
setdiff(rownames(seurat_list[[1]]), rownames(denoist_counts[[1]]))

names(denoist_counts) <- res_files

seurat_list_empty <- list()

for(xx in names(seurat_list)){
  message(xx)
  seurat_dataa <- CreateSeuratObject(counts = as.matrix(denoist_counts[[xx]]))
  
  seurat_dataa[["denoist_RNA"]] <- seurat_dataa[["RNA"]]
  DefaultAssay(seurat_dataa) <- "denoist_RNA"
  seurat_dataa[["RNA"]] <- NULL
  
  seurat_list_empty[[xx]] <- seurat_dataa
}

seurat_list_empty[[1]]
seurat_list_empty[[1]][["RNA"]]
seurat_list_empty[[1]][["denoist_RNA"]]
seurat_list_empty[[1]][["nucleus_RNA"]]
identical(rownames(seurat_list_empty[[1]][["denoist_RNA"]]), rownames(seurat_list_empty[[1]][["RNA"]]))
identical(colnames(seurat_list_empty[[1]][["denoist_RNA"]]), colnames(seurat_list_empty[[1]][["RNA"]]))

#library(scCustomize)
#merged_spatial <- Merge_Seurat_List(seurat_list)

# merged_spatial <- merge(x = seurat_list_empty[[1]], y = seurat_list_empty[2:3])
merged_spatial_denoist <- merge(x = seurat_list_empty[[1]], y = seurat_list_empty[2:length(seurat_list_empty)])

merged_spatial_denoist <- merge(x = merged_spatial_denoist[[1]], y = merged_spatial_denoist[2:length(merged_spatial_denoist)])
merged_spatial_denoist <- JoinLayers(merged_spatial_denoist)

identical(rownames(cell_seurat_data[["RNA"]]), rownames(merged_spatial_denoist[["denoist_RNA"]]))
identical(colnames(cell_seurat_data[["RNA"]]), colnames(merged_spatial_denoist[["denoist_RNA"]]))

cell_seurat_data[["denoist_RNA"]] <- merged_spatial_denoist[["denoist_RNA"]]

saveRDS(cell_seurat_data, "/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_denoIST.rds")
