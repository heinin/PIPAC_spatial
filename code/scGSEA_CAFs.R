#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 03/07/2024
# Description: scGSEA on CAF scRNAseq data
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(scGSVA)
library(Seurat)
library(gprofiler2)
library(dplyr)
library(plyr)
library(tidyverse)

#==============================================================================
# scGSVA analysis
#==============================================================================

load("/home/hnatri/PIPAC_spatial/templates.CMS.rda")

seurat_data <- readRDS("/scratch/hnatri/PIPAC/PDX_CAF_Seurat_integrated.rds")

head(templates.CMS)
genenames <- gorth(templates.CMS$symbol, source_organism = "hsapiens", target_organism = "mmusculus")

templates.CMS$gene <- mapvalues(templates.CMS$symbol,
                                from = genenames$input,
                                to = genenames$ortholog_name)

hsko <- buildAnnot(species="mouse", keytype="SYMBOL", anntype="GO")
hsko@species
hsko@anntype <- "custom"
hsko@keytype

typeof(hsko@annot)
head(hsko@annot)

annots <- templates.CMS
annots <- annots[,c("gene", "class")]
colnames(annots) <- c("GeneID", "Term")
annots$Annot <- annots$Term

hsko@annot <- as.data.frame(annots)

DefaultAssay(seurat_data) <- "RNA"
seurat_data <- JoinLayers(seurat_data)

res <- scgsva(obj=seurat_data, annot=hsko, useTerm=F)

res_df <- as.data.frame(res)

saveRDS(res_df, "/scratch/hnatri/PIPAC/scGSVA_out.rds")

q(save = "no")

res_df <- readRDS("/scratch/hnatri/PIPAC/scGSVA_out.rds")

#treg_paper_markers$Annot <- gsub(" ", ".", treg_paper_markers$Annot)

print(unique(annots$Annot))

for(i in unique(annots$Annot)){
  seurat_data@meta.data[,i] <- res_df[,i]
}

# Pseudobulk of tumor
t_cells <- c(12, 13)
myeloid <- c(10, 11)
immune <- c(17, 21, 20, 22, 24, 23)
tumor_caf <- setdiff(unique(seurat_integrated$harmony_snn_res.0.8), c(t_cells, myeloid, immune))
tumor_only <- subset(seurat_integrated, subset = harmony_snn_res.0.8 %in% tumor_caf)

VlnPlot(seurat_data,
        features = c("CMS4"),
        group.by = "Sample",
        pt.size = 0,
        ncol = 1) &
  geom_boxplot(width=0.6, fill="white") &
  NoLegend() &
  xlab("") &
  ylab("CMS4 signature") &
  ggtitle("")

g1 <- seurat_data@meta.data %>%
  filter(Sample == "CT26")

g2 <- seurat_data@meta.data %>%
  filter(Sample == "CT26_CAF")

summary(g1$CMS4)
summary(g2$CMS4)
t.test(g1$CMS4, g2$CMS4)

filename <- "/home/hnatri/PIPAC_spatial/CT26_CAF_CMS4signature.pdf"

pdf(file=filename,
    width = 4,
    height = 4)

VlnPlot(seurat_data,
        features = c("CMS4"),
        group.by = "Sample",
        pt.size = 0,
        ncol = 1) &
  #geom_boxplot(width=0.6, fill="white") &
  NoLegend() &
  xlab("") &
  ylab("CMS4 signature") &
  ggtitle("")

dev.off()

VlnPlot(seurat_data,
        features = c("CMS1", "CMS2", "CMS3", "CMS4"),
        group.by = "Sample",
        pt.size = 0,
        ncol = 2) &
  geom_boxplot(width=0.6, fill="white") &
  NoLegend() &
  xlab("")
  #ylab("CMS4 signature") &
  #ggtitle("")

VlnPlot(tumor_only,
        features = c("CMS1", "CMS2", "CMS3", "CMS4"),
        group.by = "Sample",
        pt.size = 0,
        ncol = 2) &
  geom_boxplot(width=0.6, fill="white") &
  NoLegend() &
  xlab("")

VlnPlot(tumor_only,
        features = c("Cluster1", "Cluster2", "Cluster3", "Cluster4"),
        group.by = "Sample",
        pt.size = 0,
        ncol = 2) &
  geom_boxplot(width=0.6, fill="white") &
  NoLegend() &
  xlab("")

plot(tumor_only$Cluster1, tumor_only$CMS1)
