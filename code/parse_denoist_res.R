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

cell_seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat.rds")

# DenoIST results
res_files <- list.files("/scratch/hnatri/DenoIST/TMA_outputs/")

# Adjusted counts
denoist_counts <- lapply(res_files, function(tma){
  message(tma)
  res <- readRDS(paste0("/scratch/hnatri/DenoIST/TMA_outputs/", tma))
  res <- as.data.frame(res[["adjusted_counts"]])
  
  return(res)
})

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

for(xx in names(seurat_list)){
  message(xx)
  seurat_list[[xx]][["denoist_RNA"]] <- CreateAssay5Object(data = as.matrix(denoist_counts[[xx]]))
}

merged_spatial <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
merged_spatial <- JoinLayers(merged_spatial)

saveRDS(merged_spatial, "/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat_denoIST.rds")
