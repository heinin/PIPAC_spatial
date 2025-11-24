#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 08/23/2025
# Description: Denoising cellular transcripts
#==============================================================================#

#==============================================================================
# Import packages
#==============================================================================

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
  library(DenoIST)
  library(SpatialExperiment)
  library(arrow)
  library(hexbin)
  library(flexmix)
  library(pbapply)
  library(workflowr)})

#==============================================================================
# Environment variables and helper functions
#==============================================================================

setwd("/home/hnatri/PIPAC_spatial/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

source("/home/hnatri/PIPAC_spatial/code/PIPAC_colors_themes.R")
source("/home/hnatri/PIPAC_spatial/code/plot_functions.R")

#### NEIGHBOR_OFFSET.R ----
local_offset_distance_with_background <- function(mat,
                                                  tx,
                                                  coords,
                                                  tx_x = "x",
                                                  tx_y = "y",
                                                  feature_label = "gene",
                                                  distance = 50,
                                                  nbins = 200,
                                                  cl = 1) {
  message('Calculating global background...')
  # filter by qv20
  tx <- tx[tx[['qv']] >= 20,]
  
  # Create hexagonal bins
  hex_bins <- hexbin(tx[[tx_x]], tx[[tx_y]], xbins = nbins, IDs = TRUE)
  
  x_range <- diff(range(tx[, tx_x]))
  hex_radius <- x_range / hex_bins@xbins / sqrt(3)
  
  # Calculate hexbin area
  hex_area <- (3 * sqrt(3) / 2) * hex_radius^2
  
  # Assign transcripts to bins
  tx$hexbin_id <- hex_bins@cID
  #tx[,"feature_name"] <- tx[, feature_label]
  
  # Count gene occurrences per bin
  gene_bin_counts <- tx %>%
    group_by(hexbin_id, feature_name) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  # Pivot to wide gene-by-bin matrix
  gene_bin_matrix <- gene_bin_counts %>%
    pivot_wider(names_from = hexbin_id, values_from = count, values_fill = 0)
  
  # Fix column name
  #gene_bin_matrix <- as.data.frame(gene_bin_matrix)
  #names(gene_bin_matrix)[1] <- "feature_name"
  
  # Sum duplicates and format as matrix
  gene_bin_matrix <- gene_bin_matrix %>%
    #as.tibble() %>%
    group_by(feature_name) %>%
    dplyr::summarise(across(everything(), sum), .groups = "drop") %>% 
    column_to_rownames("feature_name")
  
  # Match to `mat`
  gene_bin_matrix <- gene_bin_matrix[rownames(mat), ]
  gene_bin_matrix <- gene_bin_matrix[, sapply(gene_bin_matrix, is.numeric)]
  
  # Clean invalid values
  gene_bin_matrix[!is.finite(as.matrix(gene_bin_matrix))] <- 0
  gene_bin_matrix[is.na(gene_bin_matrix)] <- 0
  
  # Compute total per bin
  bin_total <- colSums(gene_bin_matrix)
  
  # Check variability
  if (length(unique(bin_total)) <= 1) {
    stop("GMM cannot be fit: 'bin_total' is constant or lacks sufficient variability")
  }
  
  # Fit GMM
  message("Running GMM...")
  mo1 <- FLXMRglm(family = "gaussian")
  mo2 <- FLXMRglm(family = "gaussian")
  flexfit <- tryCatch({
    flexmix(x ~ 1, data = data.frame(x = bin_total), k = 2, model = list(mo1, mo2))
  }, error = function(e) {
    stop("flexmix failed during GMM fit: ", e$message)
  })
  
  # Get cluster parameters
  c1 <- parameters(flexfit, component = 1)[[1]]
  c2 <- parameters(flexfit, component = 2)[[1]]
  gmm_means <- c(c1[1], c2[1])
  smaller_mean_component <- which.min(gmm_means)
  
  # Extract empty bins
  empty_bin_matrix <- gene_bin_matrix[, clusters(flexfit) == smaller_mean_component, drop = FALSE]
  empty_bin_matrix <- empty_bin_matrix[, colSums(empty_bin_matrix) > 0, drop = FALSE]
  
  # Compute background offset
  per_unit_sum <- rowSums(empty_bin_matrix) / (ncol(empty_bin_matrix) * hex_area)
  scaled_sum <- per_unit_sum * distance^2 * pi
  bg_offset <- ifelse(scaled_sum == 0, 1, ceiling(scaled_sum))
  
  # Neighbor-finding function
  get_neighbors_within_distance <- function(coords, distance) {
    coords_mat <- as.matrix(coords)
    mode(coords_mat) <- "numeric"
    neighbors <- pblapply(seq_len(nrow(coords)), function(i) {
      dists <- sqrt(rowSums2((coords_mat - coords_mat[i, ])^2))
      which(dists <= distance)
    }, cl = cl)
    return(neighbors)
  }
  
  message("Finding neighbours...")
  neighbors <- get_neighbors_within_distance(coords[, c(1,2)], distance)
  
  get_local_offset <- function(idx, neighbors, mat) {
    if (length(neighbors[[idx]]) == 0) {
      offset <- rep(0, nrow(mat)) + mat[, idx]
    } else if (length(neighbors[[idx]]) == 1) {
      offset <- mat[, neighbors[[idx]]] + mat[, idx]
    } else {
      offset <- rowSums2(mat[, neighbors[[idx]]]) + mat[, idx]
    }
    return(offset)
  }
  
  message("Calculating local offset...")
  res <- pblapply(seq_len(ncol(mat)), get_local_offset, neighbors, mat, cl = cl)
  res_mat <- do.call(cbind, res)
  colnames(res_mat) <- colnames(mat)
  
  # Add global background offset
  res_mat <- sweep(res_mat, 1, bg_offset, "+")
  
  return(res_mat)
}

#### PMM_MODEL.R ----
solve_poisson_mixture <- function(x, s,
                                  max_iter = 5000,
                                  tol = 1e-6,
                                  pi_inits = runif(10, min = 0, max = 0.5),
                                  verbose = FALSE) {
  
  n <- length(x)
  
  # Store indices of non-zero s
  non_zero_indices <- which(s > 0)
  
  # Remove entries with s = 0
  x <- x[non_zero_indices]
  s <- s[non_zero_indices]
  
  best_result <- NULL
  best_log_lik <- -Inf
  
  for (pi_init in pi_inits) {
    # Initialize parameters
    lambda1 <- mean(x) / mean(s)
    lambda2 <- mean(x) / (2 * mean(s))
    pi <- pi_init
    
    if (verbose) {
      cat("Initial parameters for pi =", pi_init, ":\n")
      cat("lambda1:", lambda1, "lambda2:", lambda2, "pi:", pi, "\n")
    }
    
    log_likelihood <- function(x, s, lambda1, lambda2, pi) {
      sum(log(pi * dpois(x, s * lambda1) + (1 - pi) * dpois(x, s * lambda2)))
    }
    
    log_lik <- log_likelihood(x, s, lambda1, lambda2, pi)
    
    for (iter in 1:max_iter) {
      # E-step: calculate responsibilities
      tau1 <- pi * dpois(x, s * lambda1)
      tau2 <- (1 - pi) * dpois(x, s * lambda2)
      gamma <- tau1 / (tau1 + tau2)
      
      # M-step: update parameters
      lambda1 <- sum(gamma * x) / sum(gamma * s)
      lambda2 <- sum((1 - gamma) * x) / sum((1 - gamma) * s)
      pi <- mean(gamma)
      
      if (verbose) {
        cat("Iteration", iter, "parameters:\n")
        cat("lambda1:", lambda1, "lambda2:", lambda2, "pi:", pi, "\n")
      }
      
      # Check for convergence
      new_log_lik <- log_likelihood(x, s, lambda1, lambda2, pi)
      if (!is.finite(new_log_lik) || abs(new_log_lik - log_lik) < tol) {
        if (verbose) {
          cat("Converged after", iter, "iterations\n")
          cat("Final log-likelihood:", log_lik, "\n")
        }
        break
      }
      
      log_lik <- new_log_lik
    }
    
    if (log_lik > best_log_lik) {
      best_log_lik <- log_lik
      if(abs(lambda1 - lambda2) > 1e-2) {
        # Store the best parameters
        best_result <- list(lambda1 = lambda1,
                            lambda2 = lambda2,
                            pi = pi,
                            log_lik = log_lik,
                            gamma = gamma)
      }else{
        # If model collapse occurs, keep everything
        best_result <- list(lambda1 = lambda1,
                            lambda2 = lambda2,
                            pi = pi,
                            log_lik = log_lik,
                            gamma = rep(1, length(x)))
      }
    }
  }
  
  # Assign memberships
  memberships <- ifelse(best_result$gamma > 0.6, 1, 0)
  
  # TODO: if memberships are all 0, set to 1
  if (all(memberships == 0)) {
    memberships <- rep(1, length(memberships))
  }
  
  # Pad the results to match the original input length
  full_memberships <- rep(1, n)
  full_memberships[non_zero_indices] <- memberships
  
  full_posterior <- rep(1, n)
  full_posterior[non_zero_indices] <- best_result$gamma
  
  return(list(memberships = full_memberships,
              posterior = full_posterior,
              lambda1 = best_result$lambda1,
              lambda2 = best_result$lambda2,
              pi = best_result$pi,
              log_lik = best_result$log_lik))
}

# Function to apply solve_poisson_mixture to a single column
apply_poisson_mixture_single <- function(index, c_matrix, s_matrix) {
  test_c <- c_matrix[, index]
  test_s <- s_matrix[, index]
  result <- tryCatch({
    solve_poisson_mixture(test_c, test_s)
  }, error = function(e) {
    list(memberships = rep(1, length(test_c)),
         posterior = rep(1, length(test_c)),
         lambda1 = NA, lambda2 = NA, pi = NA)
  })
  return(result)
}

#### DENOIST.R ----
utils::globalVariables(c("feature_name", "hexbin_id", "count"))

denoist <- function(mat, tx, coords = NULL,
                    tx_x = "x",
                    tx_y = "y",
                    feature_label = "gene",
                    distance = 50, nbins = 200, cl = 1, out_dir = NULL,
                    output_name){
  # TODO:check input type
  if(inherits(mat, "SpatialExperiment")){
    coords <- spatialCoords(mat)
    mat <- assay(mat)
    #remove NegControl and BLANKS
    mat <- mat[!grepl("NegControl|BLANK", rownames(mat)),]
  }else if(is.null(coords)){
    stop("coords must be provided")
  }
  
  # calculate neighbour_offset
  message("Calculating neighbour offset...")
  off_mat <- local_offset_distance_with_background(mat = mat,
                                                   coords = coords,
                                                   tx = tx,
                                                   tx_x = tx_x,
                                                   tx_y = tx_y,
                                                   feature_label = feature_label,
                                                   distance = distance,
                                                   nbins = nbins,
                                                   cl = cl)
  
  # Apply the Poisson mixture model
  message("Applying the Poisson mixture model...")
  results <- pblapply(1:ncol(mat),
                      apply_poisson_mixture_single,
                      mat,
                      off_mat,
                      cl = cl)
  
  # return neighbour_offset, adjusted_counts, posterior, params
  message("Tidying up results...")
  memberships_matrix <- do.call(cbind, lapply(results, function(res) res["memberships"]))
  memberships_matrix <- do.call(cbind, memberships_matrix)
  colnames(memberships_matrix) <- colnames(mat)
  rownames(memberships_matrix) <- rownames(mat)
  
  adjusted_counts <- mat * memberships_matrix
  colnames(adjusted_counts) <- colnames(mat)
  rownames(adjusted_counts) <- rownames(mat)
  
  # save the results
  if(!is.null(out_dir)){
    if(!dir.exists(out_dir)){
      dir.create(out_dir, recursive = TRUE)
    }
    saveRDS(list(memberships = memberships_matrix,
                 adjusted_counts = adjusted_counts,
                 params = results),
            file = paste0(out_dir, output_name))
  }
  
  return(list(memberships = memberships_matrix,
              adjusted_counts = adjusted_counts,
              params = results))
}

#==============================================================================
# Import data
#==============================================================================

# Copied to isilon /tgen_labs/banovich/PIPAC/Seurat
# /tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat.rds.rds
cell_seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_metadata.rds")

#==============================================================================
# Preparing inputs and running DenoIST
#==============================================================================

# mat : SpatialExperiment object (with the counts in assay() slot) or a count
# matrix with genes as rows and cells as columns.
# tx: Transcript data frame (a data frame with each row being an individual
# transcript, with columns specifying each transcripts' coordinates and qv).
# If your transcript file is not from Xenium and has no qv score, you can set a
# dummy column of qv = 20 for all transcripts. This workaround should not be
# needed in future updates.

# Looping through all TMAs
tma_paths <- list("MR_PIPAC-TMA1" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0034123__MR_PIPACTMA1__20240911__210055/transcripts.parquet",
                  "MR_PIPAC-TMA2" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0034221__MR_PIPACTMA2-Rerun__20241107__213721/transcripts.parquet",
                  "MR_PIPAC-TMA3" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0024823__MR_PIPACTMA3-Rerun__20241115__191438/transcripts.parquet",
                  "MR_PIPAC-TMA4" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0033982__MR_PIPACTMA4__20240920__231140/transcripts.parquet",
                  "MR_PIPAC-TMA5" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0033981__MR_PIPACTMA5__20240920__231140/transcripts.parquet",
                  "MR_PIPAC-TMA6" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0034149__MR_PIPACTMA6__20240923__231239/transcripts.parquet",
                  "MR_PIPAC-TMA7" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0034152__MR_PIPACTMA7__20240923__231239/transcripts.parquet",
                  "MR_PIPAC-TMA8" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0033448__MR_PIPACTMA8__20241224__204450/transcripts.parquet",
                  "MR_PIPAC-TMA9" = "/tgen_labs/banovich/xenium_run_folders/PIPACTMA/output-XETG00048__0041382__MR_PIPACTMA9__20250109__205856/transcripts.parquet")

# c("MR_PIPAC-TMA8", "MR_PIPAC-TMA9")
for(tma in names(tma_paths)){
  message("Running for ", tma)
  
  cell_seurat_data_tma1 <- subset(cell_seurat_data, subset = TMA == tma)
  
  # Count matrix
  mat_tma1 <- as.matrix(cell_seurat_data_tma1@assays$RNA$counts)
  
  # Coordinates
  coords_tma1 <- cell_seurat_data_tma1@meta.data %>%
    select(x_centroid, y_centroid) %>%
    rename_with(~ "x", .cols = "x_centroid") %>%
    rename_with(~ "y", .cols = "y_centroid")
  
  # Xenium data
  tx_tma1 <- read_parquet(tma_paths[[tma]])
  
  # Removing transcripts not assigned to cells
  #tx_tma1 <- tx_tma1[-which(tx_tma1$cell_id == "UNASSIGNED"),]
  
  tx_tma1_final <- tx_tma1 %>%
    rename_with(~ "x", .cols = "x_location") %>%
    rename_with(~ "y", .cols = "y_location")
  #rename_with(~ "gene", .cols = "feature_name")
  
  # Running DenoIST
  res <- denoist(mat = mat_tma1,
                 tx = tx_tma1_final,
                 coords = coords_tma1,
                 feature_label = "feature_name",
                 distance = 50, nbins = 200, cl = 1,
                 out_dir = "/scratch/hnatri/DenoIST/TMA_outputs/",
                 output_name = tma)
  
}


#res <- readRDS("/scratch/hnatri/DenoIST/TMA_outputs/denoist_results_s2_PDL017.rds")
