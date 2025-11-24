library(dplyr)
library(data.table)
library(spatstat) 
library(reshape2) 


# general params
distance_threshold = 100 # depends on dataset, units, etc 100 generally good for 3 cells 
degree_range = 45 # degrees to split 360. 30-45 usually overlaps 1 nuclei
nneighbors = 1 # n neighbors to retain in each direction

## INPUT
# X,Y COORDS: code expects columns x_centroid and y_centroid - if different change where appropriate
# sample column names 'sample', if different change where appripriate
# cell type annotation: code expects final_CT column, if different change where appripriate

#cells <- read.csv('/scratch/avannan/early_IPF_xenium/cell_proximity/spatial_prox_input_df_2025-08-07.csv') # replace data table 
#samples <- unique(cells$sample)
seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/cell_merged_spatial_filtered_splitsamples_clustered_NN30_PC50_Seurat_denoIST.rds")

seurat_data <- subset(seurat_data, subset = Arm == "Arm3")

# Tumor only
seurat_data <- subset(seurat_data, subset = Tissue == "Tumor")

# Pretreatment only
#seurat_data <- subset(seurat_data, subset = Timepoint == 0)

cells <- seurat_data@meta.data[,c("cell_id", "Sample", "Annotation", "x_centroid", "y_centroid")]
colnames(cells)[2] <- "sample"
colnames(cells)[3] <- "final_CT"
samples <- unique(cells$sample)

rm(seurat_data)

proximity_compiled <- c()
celltype_summary_compiled <- c()
for(ii in seq_along(samples)){
  sample_id <- samples[ii]
  message(sample_id)
  
  cells_sample <- subset(cells, sample == sample_id)
  celltype_summary_sample <- cells_sample %>% 
    group_by(final_CT) %>% 
    dplyr::summarise(cell_count = n())
  celltype_summary_sample$sample <- sample_id
  celltype_summary_compiled <- rbind(celltype_summary_compiled, celltype_summary_sample)
  # create sample window
  min_x <- min(cells_sample$x_centroid)
  max_x <- max(cells_sample$x_centroid)
  min_y <- min(cells_sample$y_centroid)
  max_y <- max(cells_sample$y_centroid)
  win <- owin(xrange = c(min_x, max_x), yrange = c(min_y, max_y))
  # point map
  cells_ppp <- ppp(cells_sample$x_centroid, cells_sample$y_centroid, window = win, marks = cells_sample$cell_id)
  
  # Calculate distances and angles to all other events
  dists <- sqrt(outer(cells_ppp$x, cells_ppp$x, "-")^2 + outer(cells_ppp$y, cells_ppp$y, "-")^2)
  rownames(dists) <- marks(cells_ppp)
  colnames(dists) <- marks(cells_ppp)
  dists <- reshape2::melt(dists)
  colnames(dists) <- c('anchor', 'prox', 'distance')
  
  # Calculate angles in degrees (-180 to 180)
  # degree_range = 30
  angles <- seq(0, (360 - degree_range), by = degree_range)  # 8 directions
  n_angles <- length(angles)
  
  angs <- atan2( outer(cells_ppp$y, cells_ppp$y, "-"), outer(cells_ppp$x, cells_ppp$x, "-")) * 180 / pi
  angs <- (angs + 360) %% 360
  rownames(angs) <- marks(cells_ppp)
  colnames(angs) <- marks(cells_ppp)
  
  bin_idx <- floor(angs / degree_range) %% n_angles
  rownames(bin_idx) <- marks(cells_ppp)
  colnames(bin_idx) <- marks(cells_ppp)
  
  angs <- reshape2::melt(angs)
  colnames(angs) <- c('anchor', 'prox', 'angle')
  
  bin_idx <- reshape2::melt(bin_idx)
  colnames(bin_idx) <- c('anchor', 'prox', 'angle_bin')
  
  
  proximity_table_sample <- data.frame(anchor_cell = dists$anchor, 
                                       prox_cell = dists$prox, 
                                       distance = dists$distance, 
                                       angle = angs$angle, 
                                       idx_bin = bin_idx$angle_bin)
  proximity_table_sample <- proximity_table_sample %>% 
    filter(anchor_cell != prox_cell) %>% 
    filter(distance < distance_threshold)
  
  proximity_table_sample <- proximity_table_sample %>%
    group_by(anchor_cell, idx_bin) %>%
    mutate(idx = rank(distance, ties.method = "min")) %>%
    ungroup()
  
  ## add meta data information
  proximity_table_sample$anchor_cell_ct <- cells_sample$final_CT[match(proximity_table_sample$anchor_cell, cells_sample$cell_id)]
  proximity_table_sample$prox_cell_ct <- cells_sample$final_CT[match(proximity_table_sample$prox_cell, cells_sample$cell_id)]
  proximity_table_sample$sample <-  sample_id
  
  # get only first neighbors
  proximity_table_sample <- subset(proximity_table_sample, idx == nneighbors)
  proximity_compiled <- rbind(proximity_compiled, proximity_table_sample)
}

#write.table(proximity_compiled,
#            "/home/hnatri/PIPAC_spatial/tumor_proximity_compiled.tsv",
#            quote = F, row.names = F, sep = "\t")


## OUTPUT EXPLAINED
# anchor_cell: is the cell to which neighbors are calculated for
# prox_cell: cells assigned as proximal to the anchor_cell
# distance: distance between anchor_cell and prox_cell
# angle: the angle prox_cell is in relative to anchor_cell
# idx_bin: categorical angle assignment based on angle. bins defined by degree_range
# idx: nth neighbor assignment per idx_bin


################
## ENRCIHMENT ##
################
# 1. total number of cells
total_cells <- sum(celltype_summary_compiled$cell_count)
# 2. total cell type of interest
total_celltype <- celltype_summary_compiled %>% 
  group_by(final_CT) %>% 
  summarise(n_cells_total = sum(cell_count))
colnames(total_celltype) <- c('anchor_celltype', 'n_cells_total')
# 3. total proximal cells by cell type
total_celltype_proximal <- proximity_compiled[,c('anchor_cell_ct', 'prox_cell', 'prox_cell_ct', 'sample')]
total_celltype_proximal <- total_celltype_proximal[!duplicated(total_celltype_proximal),]
total_celltype_proximal <- total_celltype_proximal %>% 
  group_by(anchor_cell_ct, prox_cell_ct) %>% 
  dplyr::summarise(n_proximal = n())
colnames(total_celltype_proximal) <- c('anchor_celltype', 'proximal_celltype', 'n_proximal')
# 4. total all proximal cells
total_proximal <- total_celltype_proximal %>% 
  group_by(anchor_celltype) %>% 
  summarise(n_proximal_total = sum(n_proximal))
# 5. combine
total_celltype_proximal$n_celltype_total <- total_celltype$n_cells_total[match(total_celltype_proximal$proximal_celltype, total_celltype$anchor_celltype)]
total_celltype_proximal$n_proximal_total <- total_proximal$n_proximal_total[match(total_celltype_proximal$anchor_celltype, total_proximal$anchor_celltype)]
total_celltype_proximal$total_cells <- total_cells

og <- total_celltype_proximal
## run fisher's exact
#total_celltype_proximal <- read.csv('/scratch/avannan/total_celltype_proximal.csv')
total_celltype_proximal$X <- NULL
celltypes <- unique(total_celltype_proximal$anchor_celltype)
enrichment_scores_compiled <- c()
celltypes <- setdiff(celltypes, c(NA))
for(celltype in celltypes){
  print(celltype)
  dt <- as.data.table(total_celltype_proximal)
  dt <- dt[anchor_celltype == celltype]
  nrow(dt)
  
  dt[, `:=`(
    a = total_cells - (n_proximal_total + n_celltype_total) + n_proximal,  # Non-overlap, non-list1, non-list2
    b = n_proximal_total - n_proximal,                    # List1 only
    c = n_celltype_total - n_proximal,                    # List2 only
    d = n_proximal                             # Overlap = proximal cells
  )]
  
  # Running Fishers test for each row
  fisher_row <- function(a, b, c, d) {
    mat <- matrix(c(a, b, c, d), nrow = 2)
    res <- fisher.test(mat, simulate.p.value = TRUE)
    c(pval = res$p.value, or = res$estimate, or05 = res$conf.int[1], or95 = res$conf.int[2])
  }
  
  dt_res <- lapply(1:nrow(dt), function(i){
    #dt[i, c("pval", "or", "or05", "or95")]
    a <- dt[i, ]$a
    b <- dt[i, ]$b
    c <- dt[i, ]$c
    d <- dt[i, ]$d
    
    as.data.table(t(fisher_row(a, b, c, d)))
  })
  dt_res <- do.call("rbind", dt_res)
  dt <- cbind(dt, dt_res)
  colnames(dt) <- gsub("or.odds ratio", "or", colnames(dt))
  
  #dt[, c("pval", "or", "or05", "or95") := as.data.table(t(fisher_row(a, b, c, d))), by = .I]
  enrich_score_ct <- dt[, .(prox_to = proximal_celltype, prox_ct = anchor_celltype, n_proximal, pval, or, or05, or95)]
  enrich_score_ct <- within(enrich_score_ct, {
    log_or = log(or)
    log_or05 = log(or05)
    log_or95 = log(or95)
  })
  # false disovery correction per celltype (can skip and apply to all)
  enrich_score_ct <- tibble::add_column(.data = enrich_score_ct, .after = 4, pvale_adjust = p.adjust(enrich_score_ct$pval, method = 'fdr'))
  enrichment_scores_compiled <- rbind(enrichment_scores_compiled, enrich_score_ct)
}

#write.table(enrichment_scores_compiled,
#            "/home/hnatri/PIPAC_spatial/enrichment_scores_compiled.tsv",
#            quote = F, row.names = F, sep = "\t")

enrichment_scores_compiled <- read.table("/home/hnatri/PIPAC_spatial/tumor_nrichment_scores_compiled.tsv",
                                         header = T)
proximity_compiled <- read.table("/home/hnatri/PIPAC_spatial/tumor_proximity_compiled.tsv",
                                 header = T)
