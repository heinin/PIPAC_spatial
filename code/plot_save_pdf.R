#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/06/2025
# Description: Saving selected plots as PDFs
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

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
  library(workflowr)
  library(patchwork)})

#==============================================================================#
# Environment variables and helper functions
#==============================================================================#

setwd("/home/hnatri/PIPAC_spatial/")
set.seed(9999)
options(future.globals.maxSize = 30000 * 1024^2)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

source("/home/hnatri/PIPAC_spatial/code/PIPAC_colors_themes.R")
source("/home/hnatri/PIPAC_spatial/code/plot_functions.R")

#==============================================================================#
# Import data
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/PIPAC_NC50_NN20_PC20_Seurat_annotated_metadata.rds")

#==============================================================================#
# Plotting
#==============================================================================#

# UMAP
umap1 <- DimPlot(seurat_data,
                 group.by = "Annotation",
                 cols = pipac_celltype_col,
                 reduction = "umap",
                 raster = T,
                 label = T,
                 label.box = T,
                 repel = T,
                 label.size = 2) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  NoLegend() +
  ggtitle("Cell type") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  manuscript_theme

umap1

filename <- "/home/hnatri/PIPAC_spatial/annotation_dimplot.pdf"
pdf(file = filename,
    width = 4,
    height = 3)

umap1

dev.off()

# Cell type proportion scatter plot
prop_scatter <- create_clusterpropplot(seurat_object = seurat_data,
                             group_var = "Tissue",
                             group1 = "Normal",
                             group2 = "Tumor",
                             plot_var = "Annotation",
                             plot_colors = pipac_celltype_col,
                             var_names = c("Normal", "Tumor"),
                             legend_title = "") +
  umap1

prop_scatter

filename <- "/home/hnatri/PIPAC_spatial/tissue_celltypeprop_scatter.pdf"
pdf(file = filename,
    width = 3,
    height = 3)

prop_scatter

dev.off()


# Patient dimplots
patient_dimplots <- readRDS("/scratch/hnatri/PIPAC/patient_dimplots.rds")

coh004 <- wrap_plots(patient_dimplots[["COH-004"]], ncol = 4)

filename <- "/home/hnatri/PIPAC_spatial/coh004.pdf"
pdf(file = filename,
    width = 10,
    height = 10)

coh004

dev.off()
