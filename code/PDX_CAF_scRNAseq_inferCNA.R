#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/04/2026
# Description: InferCNA on CAF scRNAseq data
#==============================================================================#

#==============================================================================
# Import packages
#==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(anndata)
  library(infercna)
  library(SCEVAN)})

#==============================================================================
# Environment variables and helper functions
#==============================================================================

# Environment variables and helper functions
setwd("/home/hnatri/PIPAC_spatial/")
set.seed(9999)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)

#==============================================================================
# Import data
#==============================================================================

seurat_integrated <- readRDS("/tgen_labs/banovich/PIPAC/CAF_scRNAseq/PDX_CAF_Seurat_annotated.rds")

sort(as.character(unique(seurat_integrated$Updated_Annotation)))

seurat_integrated$Updated_Annotation <- factor(seurat_integrated$Updated_Annotation,
                                               levels = sort(as.character(unique(seurat_integrated$Updated_Annotation))))

#==============================================================================
# Running SCEVAN
#==============================================================================

samples <- unique(seurat_integrated$Sample)
listCountMtx <- lapply(samples, function(xx){
  data_subset <- subset(seurat_integrated, subset = Sample == xx)
  LayerData(data_subset, assay = "RNA", layer = "counts")
  })

names(listCountMtx) <- samples 

countmx <- LayerData(seurat_integrated, assay = "RNA", layer = "counts")

scevan_res <- pipelineCNA(countmx, organism = "mouse", plotTree = F)

saveRDS(scevan_res, "/scratch/hnatri/scevan_res.rds")





