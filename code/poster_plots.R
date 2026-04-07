#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/04/2026
# Description: Poster plots
#==============================================================================#

# Packages
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
  library(workflowr)
  library(RColorBrewer)})

# Environment variables and helper functions
setwd("/home/hnatri/PIPAC_spatial/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

source("/home/hnatri/PIPAC_spatial/code/PIPAC_colors_themes.R")
source("/home/hnatri/PIPAC_spatial/code/plot_functions.R")

# Markers
main_plot_features <- c("PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "JAK1", "JAK2",
                        "STAT4", "STAT3", "TIGIT", "GZMB", "SELL", "CD19",
                        "CD68", "CD44", "MARCO", "APOE", "C1QB", "C1QBP",
                        "CD163", "SPP1", "VCAN", "MUC5AC", "NOTCH3", "RGS5",
                        "MS4A1", "PGA5", "PLVAP", "CCL21", "LYVE", "SPARCL1",
                        "FN1", "DCN", "LUM", "COL1A2", "PDPN", "ACTA2", "LMNA",
                        "EGR3", "TP53", "JUN", "FOS", "BRAF", "KIT", "SOX9",
                        "RNF43", "EPCAM", "PECAM", "CEACAM", "TGFB1", "JUNB",
                        "TCF7", "GZMA", "GZMK", "CD47", "CD74", "LYZ", "TREM2",
                        "CD247", "COL5A1", "SPARC", "CD36", "CTLA4", "RGS1",
                        "CCL2", "SELENOP", "LYVE1", "CTSD", "MMP9", "CDK1",
                        "FLT1", "RAMP2", "CCL14", "MCAM", "NRP1", "CEACAM6",
                        "IL32", "TGFBI", "VEGFA", "LAMP1", "IGFBP7", "POSTN",
                        "C1S", "C1R", "CSF1", "NOTCH1", "FOXP3", "IL2RB",
                        "TUBA1B", "S100A9", "C3", "CXCR2", "SDC1", "CD38",
                        "IGHG1", "PRDX4", "ICAM1", "GATA2", "CXCL3")

mycafs <- c("TGFB1", "ACTA2", "TAGLN") #"HOPX", "MYL9"
icafs <- c("IL6", "CXCL12", "CXCL14") #"CXCL13", 
apcafs <- c("CD74", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DRA", "HLA-DRB1")
ficafs <- c("SFRP2", "SFRP4", "COMP", "CXCL14", "PLA2G2A", "ASPN")

main_plot_features <- c(main_plot_features, mycafs, icafs, apcafs, ficafs)
main_plot_features <- unique(main_plot_features)

# Cell type markers by group
gs4_deauth()
markers_by_group  <- gs4_get("https://docs.google.com/spreadsheets/d/1w0wrL-5KwNEWi_VpmriN7YmfwpBPlEYX0MKlwI2ZQu8/edit?usp=sharing")
markers_by_group <- read_sheet(markers_by_group, sheet = "Markers by group")

metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1sXXwOreLxjMSUoPt79c6jmaQpluWkaxA5P5HfDsed3I/edit?usp=sharing")
markers <- read_sheet(metadata, sheet = "Markers")
#final_annotations <- read_sheet(metadata, sheet = "Cell type annotations, all transcripts, final")
final_annotations <- read_sheet(metadata, sheet = "CT update clean")

# Import data
seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/Freeze/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_denoIST_annotated_updated.rds")

arm2 <- subset(seurat_data, subset = Arm == "Arm2")
arm2_sample_info <- arm2@meta.data %>% dplyr::select(c("Sample", "Patient_ID", "Tissue",
                                                       "Timepoint", "TMA", "Location_Quadrant",
                                                       "GENDER", "ETHNICITY", "race", "DISEASESITE",
                                                       "DXHISTOLOGY", "No_of_PIPACs")) %>%
  distinct() %>%
  remove_rownames()

#write.table(arm2_sample_info, "/home/hnatri/arm2_sample_info.tsv", row.names = F, sep = "\t")

arm2_sample_info %>% filter(Patient_ID == "COH-004", Timepoint == 0)

arm2_sample_info %>% group_by(Patient_ID) %>%
  dplyr::summarize(n_samples = n())

sort(unique(seurat_data$Celltype_Final))

DefaultAssay(seurat_data) <- "RNA"
seurat_data <- NormalizeData(seurat_data)
seurat_data <- ScaleData(seurat_data)

DefaultAssay(seurat_data) <- "denoist_RNA"
seurat_data <- NormalizeData(seurat_data)
seurat_data <- ScaleData(seurat_data)

# Plotting annotations

