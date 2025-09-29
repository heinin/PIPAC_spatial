#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/19/2025
# Description: Running niche construction with different parameters
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

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
  library(workflowr)
  library(patchwork)
  library(RANN)})

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

# Modified BuildNicheAssay()
mod_BuildNicheAssay <- function(
    object,
    fov,
    group.by,
    assay = "niche",
    cluster.name = "niches",
    neighbors.k = 20
) {
  # initialize an empty cells x groups binary matrix
  cells <- Cells(object[[fov]])
  group.labels <- unlist(object[[group.by]][cells, ])
  groups <- sort(unique(group.labels))
  cell.type.mtx <- matrix(
    data = 0,
    nrow = length(cells),
    ncol = length(groups)
  )
  rownames(cell.type.mtx) <- cells
  colnames(cell.type.mtx) <- groups
  # populate the binary matrix 
  cells.idx <- seq_along(cells)
  group.idx <- match(group.labels, groups)
  cell.type.mtx[cbind(cells.idx, group.idx)] <- 1
  
  # find neighbors based on tissue position
  coords <- GetTissueCoordinates(object[[fov]], which = "centroids")
  rownames(coords) <- coords[["cell"]]
  coords <- as.matrix(coords[ , c("x", "y")])
  neighbors <- FindNeighbors(coords, k.param = neighbors.k, compute.SNN = FALSE)
  
  # create niche assay
  sum.mtx <- as.matrix(neighbors[["nn"]] %*% cell.type.mtx)
  niche.assay <- CreateAssayObject(counts = t(sum.mtx))
  object[[assay]] <- niche.assay
  DefaultAssay(object) <- assay
  
  return(object)
}

#==============================================================================#
# Import data
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/PIPAC/Seurat/PIPAC_NC50_NN20_PC20_Seurat_annotated_metadata_ncells3k_nneighbors20.rds")

#==============================================================================#
# Adjusting coordinates
#==============================================================================#

# Shifting coordinates for each TMA
max(seurat_data$x_centroid) #12000
max(seurat_data$y_centroid) #22000

DimPlot(seurat_data,
        split.by = "TMA",
        group.by = "Annotation",
        cols = pipac_celltype_col,
        ncol = 5,
        reduction = "sp") +
  coord_fixed()

metadata <- seurat_data@meta.data
metadata <- metadata %>%
  mutate(adj_x_centroid = ifelse(TMA %in% c("MR_PIPAC-TMA1", "MR_PIPAC-TMA6"), (-1*x_centroid),
                                 ifelse(TMA == "MR_PIPAC-TMA2", ((x_centroid+12000)*-1),
                                        ifelse(TMA == "MR_PIPAC-TMA3", ((x_centroid+(12000*2))*-1),
                                               ifelse(TMA == "MR_PIPAC-TMA4", ((x_centroid+(12000*3))*-1),
                                                      ifelse(TMA == "MR_PIPAC-TMA5", ((x_centroid+(12000*4))*-1),
                                                                    ifelse(TMA == "MR_PIPAC-TMA7", ((x_centroid+(12000))*-1),
                                                                           ifelse(TMA == "MR_PIPAC-TMA8", ((x_centroid+(12000*2))*-1),
                                                                                  ifelse(TMA == "MR_PIPAC-TMA9", ((x_centroid+(12000*3))*-1), NA)))))))),
         adj_y_centroid = ifelse(TMA %in% c("MR_PIPAC-TMA1", "MR_PIPAC-TMA2", "MR_PIPAC-TMA3", "MR_PIPAC-TMA4", "MR_PIPAC-TMA5"), y_centroid,
                                 ifelse(TMA %in%c("MR_PIPAC-TMA6", "MR_PIPAC-TMA7", "MR_PIPAC-TMA8", "MR_PIPAC-TMA9"), y_centroid+22000, NA)))

seurat_data@meta.data <- metadata

# Add spatial dimension reduction object separately
position_xy <- cbind((metadata$adj_x_centroid),
                     (metadata$adj_y_centroid))
row.names(position_xy) <- row.names(metadata)
colnames(position_xy) <- c("SP_1", "SP_2")
seurat_data[["sp_adj"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(seurat_data))

DimPlot(seurat_data,
        group.by = "TMA",
        #cols = pipac_celltype_col,
        reduction = "sp_adj") +
  coord_fixed()

#==============================================================================#
# Excluding samples with few cells
#==============================================================================#

table(seurat_data$Sample) %>%
  as.data.frame() %>%
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  RotatedAxis() +
  xlab("") +
  ylab("# cells")

# Keeping samples with at least 5000 cells
keep_samples <- table(seurat_data$Sample) %>%
  as.data.frame() %>%
  filter(Freq >= 3000) %>%
  dplyr::select(Var1) %>%
  unlist()

seurat_data <- subset(seurat_data, subset = Sample %in% keep_samples)

#==============================================================================#
# Checking nearest neighbors for each cell
#==============================================================================#

coords <- seurat_data@meta.data[, c("adj_x_centroid", "adj_y_centroid")]
rownames(coords) <- colnames(seurat_data)

# Using RANN
coords <- as.matrix(coords[ , c("adj_x_centroid", "adj_y_centroid")])

nns <- nn2(coords, coords,
           k = 20,
           treetype = "kd",
           searchtype = "standard",
           radius = 0,
           eps = 0)

cellnames <- as.data.frame(matrix(nrow = nrow(nns$nn.idx),
                                  ncol = 2))
colnames(cellnames) <- c("cellname", "index")
cellnames$cellname <- rownames(coords)
cellnames$index <- seq(1, nrow(cellnames))

rownames(nns$nn.idx) <- rownames(coords)
colnames(nns$nn.idx) <- paste0("neighbor", seq(1, 20))

neighbors_mx <- nns$nn.idx

neighbors_df <- neighbors_mx %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellname") %>%
  pivot_longer(cols = colnames(neighbors_mx),
               names_to = "neighbor",
               values_to = "nn_index")

neighbors_df <- merge(neighbors_df, cellnames, by.x = "nn_index", by.y = "index")

neighbors_df$cell1_sample <- sapply(strsplit(neighbors_df$cellname.x, "_"), `[`, 1)
neighbors_df$cell2_sample <- sapply(strsplit(neighbors_df$cellname.y, "_"), `[`, 1)

# There are some cells with neighboring cells belonging to other samples
mismatches <- neighbors_df %>%
  filter(cell1_sample != cell2_sample)

#==============================================================================#
# Niche construction
#==============================================================================#

# Constructing the niche assay separately for each sample
samples <- unique(seurat_data$Sample)

sample_niches <- lapply(samples, function(sample){
  seurat_subset <- subset(seurat_data, subset = Sample == sample)
  coords <- seurat_subset@meta.data[, c("adj_x_centroid", "adj_y_centroid")]
  seurat_subset[["fov"]]  <- CreateFOV(coords,
                                       assay = "sp_adj",
                                       type = "centroids")
  
  seurat_subset <- mod_BuildNicheAssay(object = seurat_subset,
                                       fov = "fov",
                                       group.by = "Annotation",
                                       neighbors.k = 20)
  
  DefaultAssay(seurat_subset) <- "RNA"
  
  seurat_subset[["niche"]]
})

# Placing the niche assay
merged_niches <- merge(x = sample_niches[[1]],
                       y = sample_niches[2:length(sample_niches)])

seurat_with_niches <- seurat_data
seurat_with_niches[["niche"]] <- merged_niches

# Clustering the niche assay
seurat_with_niches <- ScaleData(seurat_with_niches,
                         assay = "niche")

for(k in seq(9, 12)){
  niche_column <- paste0("ncells3k_niche_k", k, "_n20")
  kmeans_res <- kmeans(
    x = t(seurat_with_niches[["niche"]]@scale.data),
    centers = k,
    nstart = 1)
  seurat_with_niches[[niche_column]] <- kmeans_res[["cluster"]]
}

# Saving the Seurat object with niche annotations
saveRDS(seurat_with_niches, "/tgen_labs/banovich/PIPAC/Seurat/PIPAC_NC50_NN20_PC20_Seurat_annotated_metadata_ncells3k_nneighbors20.rds")

