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

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/PIPAC_NC50_NN20_PC20_Seurat_annotated_metadata_niches.rds")

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

# Cell type proportion barplots by timepoint
# seurat_object, plot_var, group_var, group_levels, plot_levels, plot_colors, var_names, legend_title
prop_barp_arm2 <- create_barplot(seurat_object = subset(seurat_data, subset = Arm == "Arm2" & Tissue == "Tumor"),
                                 group_var = "Timepoint",
                                 group_levels = c("0", "6", "12"),
                                 plot_var = "Annotation",
                                 plot_levels = sort(unique(seurat_data$Annotation)),
                                 plot_colors = pipac_celltype_col,
                                 var_names = c("%", "Timepoint"),
                                 legend_title = "")

prop_barp_arm2

filename <- "/home/hnatri/PIPAC_spatial/timepoint_celltypeprop_barplot_arm2.pdf"
pdf(file = filename,
    width = 4.5,
    height = 5)

prop_barp_arm2

dev.off()

prop_barp_arm3 <- create_barplot(seurat_object = subset(seurat_data, subset = Arm == "Arm3" & Tissue == "Tumor"),
                                 group_var = "Timepoint",
                                 group_levels = c("0", "6"),
                                 plot_var = "Annotation",
                                 plot_levels = sort(unique(seurat_data$Annotation)),
                                 plot_colors = pipac_celltype_col,
                                 var_names = c("%", "Timepoint"),
                                 legend_title = "")

prop_barp_arm3

filename <- "/home/hnatri/PIPAC_spatial/timepoint_celltypeprop_barplot_arm3.pdf"
pdf(file = filename,
    width = 4,
    height = 5)

prop_barp_arm3

dev.off()

# Niche proportions
niche_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- 5)
names(niche_col) <- seq(1, 5)

niche_barp_arm2 <- create_barplot(seurat_object = subset(seurat_data, subset = Arm == "Arm2" & Tissue == "Tumor"),
                                  group_var = "Timepoint",
                                  group_levels = c("0", "6", "12"),
                                  plot_var = "niche_k5_n30",
                                  plot_levels = sort(unique(seurat_data$niche_k5_n30)),
                                  plot_colors = niche_col,
                                  var_names = c("%", "Timepoint"),
                                  legend_title = "")

niche_barp_arm3 <- create_barplot(seurat_object = subset(seurat_data, subset = Arm == "Arm3" & Tissue == "Tumor"),
                                  group_var = "Timepoint",
                                  group_levels = c("0", "6"),
                                  plot_var = "niche_k5_n30",
                                  plot_levels = sort(unique(seurat_data$niche_k5_n30)),
                                  plot_colors = niche_col,
                                  var_names = c("%", "Timepoint"),
                                  legend_title = "")

filename <- "/home/hnatri/PIPAC_spatial/timepoint_nicheprop_barplot_arm2_arm3.pdf"
pdf(file = filename,
    width = 7,
    height = 5)

niche_barp_arm2 + niche_barp_arm3

dev.off()

# Niche composition
niche_ct_prop <- table(seurat_data$Annotation, seurat_data$niche_k5_n30)
niche_ct_prop <- prop.table(niche_ct_prop, margin = 2)
niche_ct_prop <- niche_ct_prop %>%
  as.data.frame()

# Adding lineage
# Lineage info
gs4_deauth()
metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1sXXwOreLxjMSUoPt79c6jmaQpluWkaxA5P5HfDsed3I/edit?usp=sharing")
lineage <- read_sheet(metadata, sheet = "Cell type annotations")

niche_ct_prop$lineage <- mapvalues(niche_ct_prop$Var1,
                                   from = lineage$Annotation,
                                   to = lineage$Lineage)

niche_ct_prop$Var1 <- factor(niche_ct_prop$Var1,
                             levels = rev(unique(niche_ct_prop$Var1)))

p1 <- ggplot(niche_ct_prop, aes(x = Var2, y = Var1)) +
  geom_point(aes(size = Freq, 
                 #color = value,
                 fill = Var1), 
             stroke = 0.5,
             shape = 21) +
  guides(color = FALSE) +
  scale_fill_manual(values = pipac_celltype_col) +
  theme_bw() +
  xlab("") +
  ylab("") +
  facet_grid(lineage ~ ., scales = "free", space = "free")

filename <- "/home/hnatri/PIPAC_spatial/niche_composition_k5_n30.pdf"
pdf(file = filename,
    width = 3,
    height = 8)

p1 + NoLegend()

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
