#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 25/09/2025
# Description: Pairwise proximity by sample
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

suppressPackageStartupMessages({
  library(workflowr)
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
  library(patchwork)
  library(spatstat)})

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

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat_denoIST_metadata_ncells3k_nk20_niches.rds")

#==============================================================================#
# Overall proximity
#==============================================================================#

cell_types <- unique(seurat_data$Annotation)
celltype_pairs <- combn(cell_types, 2, simplify = FALSE) 
n <- length(celltype_pairs)
pb <- txtProgressBar(min = 1, max = n, style = 3)
proximity_results <- c()

length(celltype_pairs)

for(ii in seq_along(celltype_pairs)[1:52]){
  message(celltype_pairs[[ii]])
  
  ct_pairs <- celltype_pairs[[ii]]
  ct_pairs_out <- paste(ct_pairs, collapse = "-")
  celltypeA <- ct_pairs[1]
  celltypeB <- ct_pairs[2]
  samples <- unique(seurat_data$Sample)
  
  for(sample_id in samples){
    message(sample_id)
    
    cells_sample <- subset(seurat_data, Sample == sample_id)
    # sample_affect_status <- unique(cells_sample$sample_affect)
    if(sum(cells_sample$Annotation == celltypeA) < 1 | sum(cells_sample$Annotation == celltypeB) < 1){
      message('no cells found in sample: ', sample_id)
      results_df <- data.frame(sample = sample_id,
                               pair = ct_pairs_out,
                               binary_S = NA, 
                               normalized_T = NA)
      proximity_results <- rbind(proximity_results, results_df)
      next
    }
    
    # create sample window
    min_x <- min(cells_sample$x_centroid)
    max_x <- max(cells_sample$x_centroid)
    min_y <- min(cells_sample$y_centroid)
    max_y <- max(cells_sample$y_centroid)
    win <- owin(xrange = c(min_x, max_x), yrange = c(min_y, max_y))
    # Subset by celltype pairs
    cells_anchor <- subset(cells_sample, Annotation == celltypeA)
    cells_other <- subset(cells_sample, Annotation == celltypeB)
    X_anchor <- ppp(cells_anchor$x_centroid, cells_anchor$y_centroid, window = win, marks = cells_anchor$cell_id)
    X_other <- ppp(cells_other$x_centroid, cells_other$y_centroid, window = win, marks = cells_other$cell_id)
    # Calculate distances and angles to all other events
    dists <- sqrt(outer(X_anchor$x, X_other$x, "-")^2 + outer(X_anchor$y, X_other$y, "-")^2)
    rownames(dists) <- marks(X_anchor)
    colnames(dists) <- marks(X_other)
    
    # compute summary stats binary
    dists_binary <- dists
    d <- 50 # distance threshold 
    dists_binary[dists_binary<d] = 1
    dists_binary[dists_binary>d] = 0
    P=sum(as.matrix(dists_binary))
    S=0.5 * (P/nrow(dists_binary) + P/ncol(dists_binary))
    
    # compute summary stat for none binary dist matrix
    D_ref <- median(as.matrix(dists))
    dists_normalized <- dists/D_ref
    n_A <- nrow(dists_normalized)
    n_B <- ncol(dists_normalized)
    min_A_to_B <- apply(dists_normalized, 1, min)
    avg_A_to_B <- mean(min_A_to_B)
    min_B_to_A <- apply(dists_normalized, 2, min)
    avg_B_to_A <- mean(min_B_to_A)
    T <- 0.5 * (avg_A_to_B + avg_B_to_A)
    
    results_df <- data.frame(sample = sample_id,
                             pair = ct_pairs_out,
                             binary_S = S, 
                             normalized_T = T)
    proximity_results <- rbind(proximity_results, results_df)
    # progress bar
    Sys.sleep(0.01)
    setTxtProgressBar(pb, ii)
  }
}

write.table(proximity_results, "/scratch/hnatri/PIPAC/summarize_tissue_cell-cell_interactions_d50_chunk1.tsv",
            quote = F, sep = "\t", row.names = F)

#==============================================================================#
# 
#==============================================================================#


