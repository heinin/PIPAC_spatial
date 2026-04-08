#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/04/2026
# Description: Poster plots
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

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

#==============================================================================#
# Marker info
#==============================================================================#

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

#==============================================================================#
# Loading libraries
#==============================================================================#

# Import data
seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/Freeze/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_denoIST_annotated_updated.rds")

sort(unique(seurat_data$Celltype_Final))

DefaultAssay(seurat_data) <- "RNA"
seurat_data <- NormalizeData(seurat_data)
seurat_data <- ScaleData(seurat_data)

#==============================================================================#
# Plotting annotations
#==============================================================================#

celltype_final <- as.factor(sort(unique(seurat_data$Celltype_Final)))
celltype_final_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(celltype_final))
names(celltype_final_col) <- levels(celltype_final)

DimPlot(seurat_data,
        group.by = "Celltype_Final",
        cols = celltype_final_col,
        reduction = "umap",
        raster = T,
        label = T,
        label.box = T) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  NoLegend()

unique(seurat_data$Celltype_Final)

seurat_data$lineage <- seurat_data$Celltype_Final

seurat_data$lineage[grep("Peri", seurat_data$lineage)] <- "Perithelial"
seurat_data$lineage[grep("Endo", seurat_data$lineage)] <- "Endothelial"
seurat_data$lineage[grep("Plasma", seurat_data$lineage)] <- "Plasma"
seurat_data$lineage[grep("Treg", seurat_data$lineage)] <- "Lymphoid"
seurat_data$lineage[grep("B1", seurat_data$lineage)] <- "Lymphoid"
seurat_data$lineage[grep("Unclassified", seurat_data$lineage)] <- "Unclassified"
seurat_data$lineage[grep("CellCycle", seurat_data$lineage)] <- "Cell Cycle"
seurat_data$lineage[grep("Tumor", seurat_data$lineage)] <- "Tumor"
seurat_data$lineage[grep("^F", seurat_data$lineage)] <- "Fibroblast"
seurat_data$lineage[grep("^L", seurat_data$lineage)] <- "Lymphoid"
seurat_data$lineage[grep("^M", seurat_data$lineage)] <- "Myeloid"
seurat_data$lineage[grep("Neut", seurat_data$lineage)] <- "Myeloid"

seurat_data$lineage <- factor(seurat_data$lineage,
                              levels = c("Unclassified", "Tumor", "Cell Cycle", "Fibroblast",
                                         "Endothelial", "Perithelial",
                                         "Lymphoid", "Plasma","Myeloid"))

lineage <- as.factor(levels(seurat_data$lineage))
lineage_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(lineage))
names(lineage_col) <- levels(lineage)

lineage_col[["Unclassified"]] <- "gray90"
lineage_col[["Perithelial"]] <- "purple4"

p1 <- DimPlot(seurat_data,
        group.by = "lineage",
        cols = lineage_col,
        reduction = "umap",
        raster = T) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  ggtitle("") +
  NoAxes()

filename <- "/home/hnatri/PIPAC_spatial/lineage_UMAP.pdf"
pdf(file = filename,
    width = 6,
    height = 4)
p1
dev.off()

selected_markers <- c("CCL21", "SPP1", "CD3E", "IL7R", "CD47", "JCHAIN",
                      "IGHG1", "CXCL10", "CXCL9", "ENTPD1", "NOTCH3", "PTEN",
                      "LYVE1", "S100A8", "LYZ", "VSIR", "ACTB", "CD28", "CD274",
                      "C1QC", "C1QB", "CD163", "LGMN", "EPCAM",
                      "CEACAM6", "VIM", "PDPN", "COL1A2", "COL5A1", "ACTA2",
                      "CCL21", "RGS5", "SELL", "CD8A", "CD4", "FN1",
                      "FOXP3", "TIGIT", "STAT3", "MSLN", "LRRN4", "MARCO",
                      "FAP", "POSTN", "MKI67", "CD248", "UBE2C", "HSPH1",
                      "SELENOP", "S100A9", "SLC40A1", "TGFBI", "APOE", "CD9",
                      "DDR1", "HAVCR2", "CTLA4", "TSPAN1", "CDKN1A", "CCNE2",
                      "SIRPA", "PTPN6", "SIGLEC10", "TNFRSF14", "ITGAV",
                      "VEGFA", "NRP1", "CD36", "C1R", "CD14", "CD27", "IGKC",
                      "BANK1", "NOTCH2", "JAK2", "FLT1", "IL7R", "DGKA",
                      "CXCR4", "KLRK1", "XBP1")

seurat_data <- subset(seurat_data, subset = lineage == "Unclassified", invert = T)

seurat_data$lineage <- factor(seurat_data$lineage,
                              levels = c("Tumor", "Cell Cycle", "Fibroblast",
                                         "Endothelial", "Perithelial",
                                         "Lymphoid", "Plasma","Myeloid"))

mini_features <- c("PTPRC", "CD68", "C1QC", "SPP1", "JCHAIN", "CD4",
                   "CD3E", "CD8A", 
                   "PDGFRB", "RGS5", "PDPN", 
                   "FN1", "DCN",
                   "COL5A1",
                   "CXCL12", "MKI67", "TGFBI", 
                   "EPCAM", "CEACAM6")

p2 <- DotPlot(seurat_data,
        group.by = "lineage",
        features = mini_features,
        cols = c("gray89", "tomato3")) +
  RotatedAxis()


filename <- "/home/hnatri/PIPAC_spatial/lineage_dotplot.pdf"
pdf(file = filename,
    width = 9,
    height = 4)
p2 + ylab("") + xlab("")
dev.off()

#==============================================================================#
# Plotting Arm 3 comparative analysis
#==============================================================================#

immune <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/immune_clustered_NN20_PC20_Seurat.rds")
cafs <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/fibro_clustered_NN20_PC20_Seurat.rds")

arm3 <- subset(seurat_data, subset = Arm == "Arm3")
arm3_immune <- subset(immune, subset = Arm == "Arm3")
arm3_cafs <- subset(cafs, subset = Arm == "Arm3")

arm3_sample_info <- arm3@meta.data %>% dplyr::select(c("Sample", "Patient_ID", "Tissue",
                                                       "Timepoint", "TMA", "Location_Quadrant",
                                                       "GENDER", "ETHNICITY", "race", "DISEASESITE",
                                                       "DXHISTOLOGY", "No_of_PIPACs")) %>%
  distinct() %>%
  remove_rownames()

rm(seurat_data)

arm3_sample_info %>% group_by(Patient_ID) %>%
  dplyr::summarize(n_samples = n())

DefaultAssay(arm3) <- "RNA"
DefaultAssay(arm3_immune) <- "RNA"
DefaultAssay(arm3_cafs) <- "RNA"

arm3 <- arm3 %>% NormalizeData() %>% ScaleData()
arm3_immune <- arm3_immune %>% NormalizeData() %>% ScaleData()
arm3_cafs <- arm3_cafs %>% NormalizeData() %>% ScaleData()

#====================================
# Cell type proportions
#====================================

# For proportion testing, selecting celltypes with >500 cells in Arm 3
#keeps <- table(arm3$Celltype_Final) %>%
#  as.data.frame() %>%
#  filter(Freq > 500)

#arm3_subset <- subset(arm3, subset = Celltype_Final %in% keeps$Var1)
arm3_normal <- subset(arm3, subset = Tissue == "Normal")
arm3_tumor <- subset(arm3, subset = Tissue == "Tumor")

#arm3_immune <- subset(arm3_immune, subset = Celltype_Final %in% keeps$Var1)
#arm3_cafs <- subset(arm3_cafs, subset = Celltype_Final %in% keeps$Var1)

# Tumor immune
arm3_tumor_immune <- subset(arm3_immune, subset = Tissue == "Tumor")
keeps <- table(arm3_tumor_immune$Celltype_Final) %>%
  as.data.frame() %>%
  filter(Freq > 300)

arm3_tumor_immune <- subset(arm3_tumor_immune, subset = Celltype_Final %in% keeps$Var1)

# Normal immune
arm3_normal_immune <- subset(arm3_immune, subset = Tissue == "Normal")
keeps <- table(arm3_normal_immune$Celltype_Final) %>%
  as.data.frame() %>%
  filter(Freq > 300)

arm3_normal_immune <- subset(arm3_normal_immune, subset = Celltype_Final %in% keeps$Var1)

# Tumor CAFs
arm3_tumor_caf <- subset(arm3_cafs, subset = Tissue == "Tumor")
keeps <- table(arm3_tumor_caf$Celltype_Final) %>%
  as.data.frame() %>%
  filter(Freq > 300)

arm3_tumor_caf <- subset(arm3_tumor_caf, subset = Celltype_Final %in% keeps$Var1)

# Normal CAFs
arm3_normal_caf <- subset(arm3_cafs, subset = Tissue == "Normal")
keeps <- table(arm3_normal_caf$Celltype_Final) %>%
  as.data.frame() %>%
  filter(Freq > 300)

arm3_normal_caf <- subset(arm3_normal_caf, subset = Celltype_Final %in% keeps$Var1)

## All tumors
### Cell type proportions by sample and sample type
create_barplot(arm3_subset,
               group_var = "Sample",
               plot_var = "Celltype_Final",
               plot_levels = sort((unique(arm3_subset$Celltype_Final))),
               group_levels = sort(unique(arm3_subset$Sample)),
               plot_colors = celltype_final_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Celltype")

### Lineage proportions by timepoint in tumors
p3 <- create_barplot(arm3_tumor,
               group_var = "Timepoint",
               plot_var = "lineage",
               plot_levels = sort((unique(arm3_tumor$lineage))),
               group_levels = sort(unique(arm3_tumor$Timepoint)),
               plot_colors = lineage_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Lineage")

filename <- "/home/hnatri/PIPAC_spatial/lineage_timepoint_barplot_tumor.pdf"
pdf(file = filename,
    width = 3,
    height = 3)
p3
dev.off()

### Linege proportions by timepoint in tumor-adjacent samples
create_barplot(arm3_normal,
               group_var = "Timepoint",
               plot_var = "Celltype_Final",
               plot_levels = sort((unique(arm3_normal$Celltype_Final))),
               group_levels = sort(unique(arm3_normal$Timepoint)),
               plot_colors = celltype_final_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Celltype")

p4 <- create_barplot(arm3_normal,
               group_var = "Timepoint",
               plot_var = "lineage",
               plot_levels = sort((unique(arm3_normal$lineage))),
               group_levels = sort(unique(arm3_normal$Timepoint)),
               plot_colors = lineage_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Lineage")

filename <- "/home/hnatri/PIPAC_spatial/lineage_timepoint_barplot_normal.pdf"
pdf(file = filename,
    width = 3,
    height = 3)
p4
dev.off()

#### Testing for significance with scProportionTest
prop_test <- sc_utils(arm3_tumor)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Celltype_Final",
  sample_1 = "6", sample_2 = "0",
  sample_identity = "Timepoint")

permutation_plot(prop_test) +
  ggtitle("Tumor")

prop_test <- sc_utils(arm3_normal)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Celltype_Final",
  sample_1 = "6", sample_2 = "0",
  sample_identity = "Timepoint")

permutation_plot(prop_test) +
  ggtitle("Normal")

prop_test <- sc_utils(arm3_normal_immune)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Celltype_Final",
  sample_1 = "6", sample_2 = "0",
  sample_identity = "Timepoint")

permutation_plot(prop_test) +
  ggtitle("Immune, normal")

filename <- "/home/hnatri/PIPAC_spatial/immune_timepoint_scPropTest_normal.pdf"
pdf(file = filename,
    width = 4.5,
    height = 4)
permutation_plot(prop_test) +
  ggtitle("Immune, normal")
dev.off()

prop_test <- sc_utils(arm3_tumor_immune)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Celltype_Final",
  sample_1 = "6", sample_2 = "0",
  sample_identity = "Timepoint")

permutation_plot(prop_test) +
  ggtitle("Immune, tumor")

filename <- "/home/hnatri/PIPAC_spatial/immune_timepoint_scPropTest_tumor.pdf"
pdf(file = filename,
    width = 4.6,
    height = 4)
permutation_plot(prop_test) +
  ggtitle("Immune, tumor")
dev.off()

prop_test <- sc_utils(arm3_normal_caf)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Celltype_Final",
  sample_1 = "6", sample_2 = "0",
  sample_identity = "Timepoint")

permutation_plot(prop_test) +
  ggtitle("CAF, normal")

filename <- "/home/hnatri/PIPAC_spatial/caf_timepoint_scPropTest_normal.pdf"
pdf(file = filename,
    width = 4.5,
    height = 4)
permutation_plot(prop_test) +
  ggtitle("CAF, normal")
dev.off()

prop_test <- sc_utils(arm3_tumor_caf)
prop_test <- permutation_test(
  prop_test, cluster_identity = "Celltype_Final",
  sample_1 = "6", sample_2 = "0",
  sample_identity = "Timepoint")

permutation_plot(prop_test) +
  ggtitle("CAF, tumor")

filename <- "/home/hnatri/PIPAC_spatial/caf_timepoint_scPropTest_tumor.pdf"
pdf(file = filename,
    width = 4.5,
    height = 4)
permutation_plot(prop_test) +
  ggtitle("CAF, tumor")
dev.off()

# Celltype level for immune and CAF subsets
# Tumor immune
immune_celltype_final <- as.factor(sort(unique(immune$Celltype_Final)))
immune_celltype_final_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(immune_celltype_final))
names(immune_celltype_final_col) <- levels(immune_celltype_final)

p5 <- create_barplot(arm3_tumor_immune,
                     group_var = "Timepoint",
                     plot_var = "Celltype_Final",
                     plot_levels = sort((unique(arm3_tumor_immune$Celltype_Final))),
                     group_levels = sort(unique(arm3_tumor_immune$Timepoint)),
                     plot_colors = immune_celltype_final_col,
                     var_names =  c("Frequency (%)", ""),
                     legend_title = "Lineage")

filename <- "/home/hnatri/PIPAC_spatial/immune_timepoint_barplot_tumor.pdf"
pdf(file = filename,
    width = 3,
    height = 3)
p5
dev.off()

p6 <- create_barplot(arm3_normal_immune,
                     group_var = "Timepoint",
                     plot_var = "Celltype_Final",
                     plot_levels = sort((unique(arm3_normal_immune$Celltype_Final))),
                     group_levels = sort(unique(arm3_normal_immune$Timepoint)),
                     plot_colors = immune_celltype_final_col,
                     var_names =  c("Frequency (%)", ""),
                     legend_title = "Lineage")

filename <- "/home/hnatri/PIPAC_spatial/immune_timepoint_barplot_normal.pdf"
pdf(file = filename,
    width = 3,
    height = 3)
p6
dev.off()

# CAFs
caf_celltype_final <- as.factor(sort(unique(cafs$Celltype_Final)))
caf_celltype_final_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(caf_celltype_final))
names(caf_celltype_final_col) <- levels(caf_celltype_final)

p5 <- create_barplot(arm3_tumor_caf,
                     group_var = "Timepoint",
                     plot_var = "Celltype_Final",
                     plot_levels = sort((unique(arm3_tumor_caf$Celltype_Final))),
                     group_levels = sort(unique(arm3_tumor_caf$Timepoint)),
                     plot_colors = caf_celltype_final_col,
                     var_names =  c("Frequency (%)", ""),
                     legend_title = "Lineage")

filename <- "/home/hnatri/PIPAC_spatial/cafs_timepoint_barplot_tumor.pdf"
pdf(file = filename,
    width = 3,
    height = 3)
p5
dev.off()

p6 <- create_barplot(arm3_normal_caf,
                     group_var = "Timepoint",
                     plot_var = "Celltype_Final",
                     plot_levels = sort((unique(arm3_normal_caf$Celltype_Final))),
                     group_levels = sort(unique(arm3_normal_caf$Timepoint)),
                     plot_colors = immune_celltype_final_col,
                     var_names =  c("Frequency (%)", ""),
                     legend_title = "Lineage")

filename <- "/home/hnatri/PIPAC_spatial/cafs_timepoint_barplot_normal.pdf"
pdf(file = filename,
    width = 3,
    height = 3)
p6
dev.off()

