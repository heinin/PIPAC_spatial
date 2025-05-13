#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/24/2025
# Description: Updating metadata
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(googlesheets4)
library(tidyverse)

#==============================================================================#
# Import data
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/PIPAC_NC50_NN20_PC20_Seurat_annotated_metadata_niches.rds")

gs4_deauth()
metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1sXXwOreLxjMSUoPt79c6jmaQpluWkaxA5P5HfDsed3I/edit?usp=sharing")
arm3_best_response <- read_sheet(metadata, sheet = "Arm 3 best response")
arm2_best_response <- read_sheet(metadata, sheet = "Arm 2 CEA")

arm3_best_response <- arm3_best_response %>%
  dplyr::select(c("RPN", "respgrp"))
arm2_best_response <- arm2_best_response %>%
  dplyr::select(c("RPN", "respgrp"))

best_responses <- rbind(arm3_best_response, arm2_best_response)

setdiff(best_responses$RPN, seurat_data$RPN)
setdiff(seurat_data$RPN, best_responses$RPN)

best_responses$RPN <- sprintf("%03d", best_responses$RPN)

#==============================================================================#
# Update metadata
#==============================================================================#

seurat_data$BestResponse <- mapvalues(x = seurat_data$RPN,
                                      from = best_responses$RPN,
                                      to = best_responses$respgrp)

seurat_data$BestResponse <- gsub("NA", NA, seurat_data$BestResponse)

# Saving
saveRDS(seurat_data, "/tgen_labs/banovich/PIPAC/Seurat/PIPAC_NC50_NN20_PC20_Seurat_annotated_metadata_niches.rds")
