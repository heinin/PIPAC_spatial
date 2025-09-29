#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/05/2025
# Description: PIPAC Xenium colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)

#==============================================================================
# Colors and themes
#==============================================================================

# ggplot theme
plot_theme <- theme(text = element_text(size = 10),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 10),  
                    axis.title.x = element_text(size = 10),
                    axis.title.y = element_text(size = 10))

# ggplot theme for manuscript plots
manuscript_theme <- theme(text = element_text(size = 6),
                          axis.text.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),  
                          axis.title.x = element_text(size = 6),
                          axis.title.y = element_text(size = 6))

# Colors for plotting
# Define colors for each level of categorical variables

# Cluster colors
data_clusters <- as.factor(c(0, seq(1:14)))

pipac_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(data_clusters))
names(pipac_cluster_col) <- levels(data_clusters)

# Cell-based object has 20 clusters
cell_data_clusters <- as.factor(c(0, seq(1:19)))
pipac_cluster_20_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(cell_data_clusters))
names(pipac_cluster_20_col) <- levels(cell_data_clusters)

immune_clusters <- as.factor(c(0, seq(1:21)))

pipac_immune_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(immune_clusters))
names(pipac_immune_cluster_col) <- levels(immune_clusters)

nonimmune_clusters <- as.factor(c(0, seq(1:11)))

pipac_nonimmune_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(nonimmune_clusters))
names(pipac_nonimmune_cluster_col) <- levels(nonimmune_clusters)

# Cell type colors
gs4_deauth()
metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1sXXwOreLxjMSUoPt79c6jmaQpluWkaxA5P5HfDsed3I/edit?usp=sharing")
pipac_celltypes <- read_sheet(metadata, sheet = "Cell type colors")

pipac_celltype_col <- pipac_celltypes$Color
names(pipac_celltype_col) <- pipac_celltypes$Annotation

# Tissue
tissue_col <- c("Tumor" = "brown3",
                "Normal" = "cyan4")

# Cluster annotations
#gs4_deauth()
#cluster_annot  <- gs4_get("https://docs.google.com/spreadsheets/d/127J6C4KF7uBGKUnrPuC1mcsb_wNCN6k1zXKSCbJ6q0M/edit?usp=sharing")
#cluster_annot <- read_sheet(cluster_annot, sheet = "Cluster annotation")

#carspp1_celltypes <- cluster_annot$annot
#carspp1_celltype_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(carspp1_celltypes))
#names(carspp1_celltype_col) <- carspp1_celltypes

#myeloid_clusters <- as.factor(c(0, seq(1:16)))
#myeloid_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(myeloid_clusters))
#names(myeloid_cluster_col) <- myeloid_clusters
#
#myeloid_celltypes <- c("M1_suppressive", "M6", "M9", "M2_suppressive", "M4",
#                       "M5_suppressive", "M7", "M3_suppressive",
#                       "M8_suppressive_G2MS", "M10", "M11")
#myeloid_colors <- tinter("aquamarine3", steps = 8, crop = 3,
#                         direction = "both", adjust = 0)
#myeloid_colors <- colorRampPalette(c("aquamarine1", "darkgreen"))(nb.cols <- length(myeloid_celltypes))
#names(myeloid_colors) <- myeloid_celltypes
#myeloid_colors <- myeloid_colors[myeloid_celltypes]

# Sample type
#sample_type_col <- c("Mock" = "azure3",
#                     "Pre_CAR" = "lightskyblue3",
#                     "5050_CAR" = "darkslategray2",
#                     "CAR" = "turquoise2")


