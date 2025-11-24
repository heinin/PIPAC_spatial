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

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat.rds")

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

#==============================================================================#
# Update annotations
#==============================================================================#



#==============================================================================#
# Update metadata
#==============================================================================#

seurat_data$BestResponse <- mapvalues(x = seurat_data$RPN,
                                      from = best_responses$RPN,
                                      to = best_responses$respgrp)

seurat_data$BestResponse <- gsub("NA", NA, seurat_data$BestResponse)

seurat_data$prgs_avg1 <- mapvalues(x = seurat_data$RPN,
                                      from = prgs_scores$Subject,
                                      to = as.numeric(unlist(prgs_scores$prgs_avg1)))

seurat_data$prgs_avg2 <- mapvalues(x = seurat_data$RPN,
                                   from = prgs_scores$Subject,
                                   to = as.numeric(unlist(prgs_scores$prgs_avg2)))

seurat_data$prgs_avg3 <- mapvalues(x = seurat_data$RPN,
                                   from = prgs_scores$Subject,
                                   to = as.numeric(unlist(prgs_scores$prgs_avg3)))


# Saving
saveRDS(seurat_data, "/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat_denoIST_metadata.rds")
