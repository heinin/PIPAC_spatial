#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/04/2026
# Description: Arm 2 sample summary
#==============================================================================#

# Packages
suppressPackageStartupMessages({
  library(googlesheets4)
  library(dplyr)
  library(Seurat)
  library(tidyverse)})

# Import data
seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/Freeze/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_denoIST_annotated_updated.rds")

# Adding response
gs4_deauth()
metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1sXXwOreLxjMSUoPt79c6jmaQpluWkaxA5P5HfDsed3I/edit?usp=sharing")
arm3_best_response <- read_sheet(metadata, sheet = "Arm 3 best response")
arm2_best_response <- read_sheet(metadata, sheet = "Arm 2 CEA")
prgs_scores <- read_sheet(metadata, sheet = "PRGS")

arm3_best_response <- arm3_best_response %>%
  dplyr::select(c("RPN", "respgrp"))
arm2_best_response <- arm2_best_response %>%
  dplyr::select(c("RPN", "respgrp"))

best_responses <- rbind(arm3_best_response, arm2_best_response)

setdiff(best_responses$RPN, seurat_data$RPN)
setdiff(seurat_data$RPN, best_responses$RPN)

best_responses$RPN <- sprintf("%03d", best_responses$RPN)
prgs_scores$Subject <- sprintf("%03d", prgs_scores$Subject)

seurat_data$BestResponse <- mapvalues(x = seurat_data$RPN,
                                      from = best_responses$RPN,
                                      to = best_responses$respgrp)

seurat_data$BestResponse <- gsub("NA", NA, seurat_data$BestResponse)

# Subset Arm 2
arm2 <- subset(seurat_data, subset = Arm == "Arm2")

# Sample summary
arm2_sample_info <- arm2@meta.data %>% dplyr::select(c("Sample", "Patient_ID", "Tissue",
                                                       "Timepoint", "TMA", "Location_Quadrant",
                                                       "GENDER", "ETHNICITY", "race", "DISEASESITE",
                                                       "DXHISTOLOGY", "No_of_PIPACs", "osmonths",
                                                       "pfsmonths")) %>%
  distinct() %>%
  remove_rownames()

#write.table(arm2_sample_info, "/home/hnatri/arm2_sample_info.tsv", row.names = F, sep = "\t")

# Numbers of samples per patient
arm2_sample_info %>% group_by(Patient_ID) %>%
  dplyr::summarize(n_samples = n())

arm2_sample_info %>% filter(Patient_ID == "COH-004", Timepoint == 0)

length(unique(arm2_sample_info$Patient_ID))



