#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 10/28/2025
# Description: Plotting informative cluster markers
#==============================================================================#

#==============================================================================
# Import packages
#==============================================================================

suppressPackageStartupMessages({
  library(workflowr)
  library(arrow)
  library(dplyr)
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
# Environment and helper functions
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

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_denoIST.rds")
DefaultAssay(seurat_data) <- "denoist_RNA"
seurat_data <- NormalizeData(seurat_data)

seurat_data_immune <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_immune_clustered_NN30_PC20_Seurat.rds")
DefaultAssay(seurat_data_immune) <- "denoist_RNA"
seurat_data_immune <- NormalizeData(seurat_data_immune)

seurat_data_nonimmune <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_nonimmune_clustered_NN30_PC20_Seurat.rds")
DefaultAssay(seurat_data_nonimmune) <- "denoist_RNA"
seurat_data_nonimmune <- NormalizeData(seurat_data_nonimmune)

#==============================================================================
# Plot features
#==============================================================================

lymph <- c("PTPRC", "CD3E", "CD4", "CD8A", "TIGIT", "FOXP3", "CD19", "SELL")
myel <- c("LYZ", "MARCO", "CD68", "APOE", "C1QB", "SPP1", "PTPN2")
meso <- c("MSLN", "LRRN4")
fibro <- c("FN1", "DCN", "LUM", "FAP", "PDPN", "COL1A2", "COL5A1")
tumor_epi <- c("TP53", "EPCAM")
misc <- c("CDKN1A", "CCNE2", "MKI67", "ANXA1", "CEACAM6", "RAMP2", "NOTCH2",
          "NOTCH3", "VEGFA", "CDON")

extra <- c("NKG7", "IL2RB", "KLRK1", "CD44", "IL7R", "LGMN", "SELENOP", "CD36",
           "CD74", "S100A9", "IGKC", "SPARCL1", "MZB1", "MS4A1", "LAMP1",
           "CXCL14", "JCHAIN", "TC2N", "CD40LG", "CTSD", "TREM2")

all_features <- c(lymph,
                  myel,
                  meso,
                  fibro,
                  tumor_epi,
                  misc,
                  extra)

DotPlot(seurat_data,
        group.by = "leiden_0.5",
        features = all_features,
        cols = c("gray89", "tomato3")) +
  RotatedAxis()

seurat_data_immune$leiden_1.5 <- factor(seurat_data_immune$leiden_1.5,
                                 levels = sort(unique(seurat_data_immune$leiden_1.5)))
Idents(seurat_data_immune) <- seurat_data_immune$leiden_1.5

DotPlot(seurat_data_immune,
        group.by = "leiden_1.5",
        features = all_features,
        cols = c("gray89", "tomato3")) +
  RotatedAxis()

create_dotplot_heatmap(seurat_object = seurat_data_immune,
                       plot_features = all_features,
                       group_var = "leiden_1.5",
                       group_colors = pipac_cell_immune_cluster_col,
                       column_title = "",
                       row_km = 8,
                       col_km = 8,
                       row.order = NULL,
                       col.order = NULL)

seurat_data_nonimmune$leiden_1.5 <- factor(seurat_data_nonimmune$leiden_1.5,
                                        levels = sort(unique(seurat_data_nonimmune$leiden_1.5)))
Idents(seurat_data_nonimmune) <- seurat_data_nonimmune$leiden_1.5

DotPlot(seurat_data_nonimmune,
        group.by = "leiden_1.5",
        features = all_features,
        cols = c("gray89", "tomato3")) +
  RotatedAxis()

create_dotplot_heatmap(seurat_object = seurat_data_nonimmune,
                       plot_features = all_features,
                       group_var = "leiden_1.5",
                       group_colors = pipac_cell_nonimmune_cluster_col,
                       column_title = "",
                       row_km = 8,
                       col_km = 8,
                       row.order = NULL,
                       col.order = NULL)
