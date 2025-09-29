#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/06/2025
# Description: Plotting PIPAC patient metadata
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

suppressMessages({library(Seurat)
                  library(ggplot2)
                  library(data.table)
                  library(dplyr)
                  library(patchwork)
                  library(cowplot)
                  library(tidyr)
                  library(ggrepel)
                  library(tidyverse)
                  library(googlesheets4)
                  library(ComplexHeatmap)})

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

unique(seurat_data$osmonths)

#==============================================================================#
# Plotting
#==============================================================================#

keep_cols <- c("Patient_ID",
               "age",
               "ETHNICITY",
               "race",
               "GENDER",
               "osmonths",
               "pfsmonths",
               "Arm",
               "DISEASESITE",
               "DXHISTOLOGY",
               "No_of_PIPACs",
               "Arm")

plot_data <- seurat_data@meta.data %>% dplyr::select(all_of(keep_cols)) %>%
  distinct() %>%
  arrange(Patient_ID)

# Adding the number of samples
n_samples <- seurat_data@meta.data %>% dplyr::select(all_of(c("Patient_ID", "Sample"))) %>%
  distinct() %>%
  group_by(Patient_ID) %>%
  dplyr::summarize(n_samples = n())

plot_data <- merge(plot_data, n_samples, by = "Patient_ID")
colnames(plot_data) <- gsub("age", "Age", colnames(plot_data))
colnames(plot_data) <- gsub("ETHNICITY", "Ethnicity", colnames(plot_data))
colnames(plot_data) <- gsub("race", "Race", colnames(plot_data))
colnames(plot_data) <- gsub("GENDER", "Gender", colnames(plot_data))
colnames(plot_data) <- gsub("osmonths", "Overall Survival (mo)", colnames(plot_data))
colnames(plot_data) <- gsub("pfsmonths", "Progress-free Survival (mo)", colnames(plot_data))
colnames(plot_data) <- gsub("DISEASESITE", "Disease Site", colnames(plot_data))
colnames(plot_data) <- gsub("DXHISTOLOGY", "Histological Diagnosis", colnames(plot_data))
colnames(plot_data) <- gsub("No_of_PIPACs", "Number of PIPACs", colnames(plot_data))
colnames(plot_data) <- gsub("n_samples", "Number of Samples", colnames(plot_data))

unique(plot_data$Ethnicity)
unique(plot_data$Race)
unique(plot_data$Gender)
unique(plot_data$`Disease Site`)
unique(plot_data$`Histological Diagnosis`)
unique(plot_data$`Number of PIPACs`)
unique(plot_data$Arm)

plot_data$`Overall Survival (mo)` <- unlist(plot_data$`Overall Survival (mo)`)
plot_data$`Progress-free Survival (mo)` <- unlist(plot_data$`Progress-free Survival (mo)`)

#plot_data$`Number of PIPACs` <- as.numeric(plot_data$`Number of PIPACs`)
plot_data$Age <- as.numeric(plot_data$Age)

# ComplexHeatmap
gender_col <- c("Male" = "wheat3",
                "Female" = "palegreen3")
race_col <- c("Caucasian" = "aliceblue",
              "Pac-Isl" = "antiquewhite",
              "Asian" = "aquamarine3",
              "Not disclosed" = "gray80")
ethnicity_col <- c("Non-Hispanic or Non-Latino" = "lightpink1",
                   "Hispanic or Latino" = "lemonchiffon2")
site_col <- c("Colorectal" = "darkblue",
              "Appendiceal" = "lightblue3")
arm_col <- c("Arm2" = "purple1",
             "Arm3" = "steelblue")

age_col <- colorRamp2(c(0, max(as.numeric(plot_data$Age))), c("white", "brown3"))

n_PIPACs_col <- colorRampPalette(brewer.pal(9, "Reds"))(nb.cols <- 3)
names(n_PIPACs_col) <- c("1", "2", "3")

diagnosis_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(unique(plot_data$`Histological Diagnosis`)))
names(diagnosis_col) <- sort(unique(plot_data$`Histological Diagnosis`))

survival_col <- colorRamp2(c(0, 50), c("white", "deeppink3"))

samples_col <- colorRamp2(c(0, 16), c("white", "forestgreen"))

col <- list(Gender = gender_col,
            Race = race_col,
            Ethnicity = ethnicity_col,
            "Disease Site" = site_col,
            Age = age_col,
            "Histological Diagnosis" = diagnosis_col,
            "Number of PIPACs" = n_PIPACs_col,
            os = survival_col,
            pfs = survival_col,
            "Number of Samples" = samples_col,
            "Arm" = arm_col)

ha <- HeatmapAnnotation(Patient = anno_text(as.character(plot_data$Patient_ID)),
                        Arm = plot_data$Arm,
                        Gender = plot_data$Gender,
                        Race = plot_data$Race,
                        Ethnicity = plot_data$Ethnicity,
                        Age = plot_data$Age,
                        "Disease Site" = plot_data$`Disease Site`,
                        #"Histological Diagnosis" = plot_data$`Histological Diagnosis`,
                        "Number of PIPACs" = plot_data$`Number of PIPACs`,
                        "Number of Samples" = plot_data$`Number of Samples`,
                        os = plot_data$`Overall Survival (mo)`,
                        pfs = plot_data$`Progress-free Survival (mo)`,
                        col = col,
                        na_col = "gray80",
                        annotation_name_side = "left",
                        gp = gpar(fontsize = 10, col = "white", lwd = 2),
                        annotation_legend_param = list(
                          os = list(
                            title = "Overall Survival (mo)",
                            at = c(0, 10, 20, 30, 40, 50),
                            labels = c(0, 10, 20, 30, 40, 50)),
                          pfs = list(
                            title = "Progress-free Survival (mo)",
                            at = c(0, 10, 20, 30, 40, 50),
                            labels = c(0, 10, 20, 30, 40, 50))))

mat <- matrix(nrow = 0, ncol = nrow(plot_data))
heatmap <- Heatmap(mat,
                   top_annotation = ha,
                   width = ncol(mat)*unit(5, "mm"), 
                   height = nrow(mat)*unit(5, "mm"))

draw(heatmap,
     annotation_legend_side = "bottom")

dev.off()

filename <- "/home/hnatri/PIPAC_spatial/demographics_grid.pdf"
pdf(file = filename,
    width = 8,
    height = 8)
draw(heatmap,
     annotation_legend_side = "bottom")
dev.off()

