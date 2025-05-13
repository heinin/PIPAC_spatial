mcgee_clustered_obj <- readRDS("/scratch/aoill/projects/McGee_mouse_spatial/clustered_obj_all_final_annotations_2025_01_16.rds")

# DGE analysis with potential contamination genes removed ----
## 6 weeks ----
# Subset weeks
clustered_obj_6_weeks <- subset(mcgee_clustered_obj, subset = Weeks == "6")
cell_expr <- clustered_obj_6_weeks[["RNA"]]$counts
dim(cell_expr)
cell_object_meta <- clustered_obj_6_weeks@meta.data


sample_celltype_levels <-  paste0(as.character(cell_object_meta$ID),"_",
                                  as.character(cell_object_meta$ct_final))

cell_RNA_expr_perCT <- t(sapply(by(t(cell_expr),sample_celltype_levels,colSums),identity))
dim(cell_RNA_expr_perCT)

celltype_levels <- cell_object_meta$ct_final
cell_expr_count_per1K <- apply(cell_expr,2,function(x){
  res = x/sum(x)
  res * 1000
})

prop_thresh <- t(sapply(by(t(cell_expr_count_per1K),celltype_levels, 
                           function(x){
                             colSums(x >= 3)/nrow(x) 
                           }),identity))

test_gene_prop_thresh <- data.frame(prop_thresh,
                                    check.names = FALSE) %>% 
  mutate(ct_final = rownames(prop_thresh)) %>%
  tidyr::pivot_longer(1:350,
                      names_to = "gene",
                      values_to = "prop_cell_expressing_min3")

test_gene_prop_thresh_filter0.3 <- test_gene_prop_thresh %>% 
  filter(prop_cell_expressing_min3 >=0.3)
min(table(test_gene_prop_thresh_filter0.3$ct_final))
max(table(test_gene_prop_thresh_filter0.3$ct_final))


# Loop through each cell type and perform DGE analysis on all cell types
cts_to_analyze <- levels(as.factor(clustered_obj_6_weeks@meta.data$ct_final))

all_results <- c()
all_results_efit <- c()

for (ctn in 1:length(cts_to_analyze)) {
  ct <- cts_to_analyze[ctn]
  print(ct)
  
  # get the genes to test
  test_genes <- test_gene_prop_thresh_filter0.3 %>%
    filter(ct_final == ct)
  test_genes <- test_genes$gene
  
  message(length(test_genes))
  
  
  # subset object to cell type
  obj_CT <- subset(clustered_obj_6_weeks, subset = ct_final == ct)
  
  # Extract the normalized expression data (data is log-transformed counts)
  expr_matrix <- as.data.frame(as.matrix(obj_CT@assays$RNA@data))
  expr_matrix_t <- t(expr_matrix)
  expr_matrix_t <- as.data.frame(expr_matrix_t)
  
  #all_genes <- colnames(expr_matrix_t)
  
  # Add sample information to the data
  expr_matrix_t$ID <- obj_CT@meta.data$ID
  
  # Summarize by Sample to calculate the mean expression for each gene and sample
  mean_expression_per_sample <- expr_matrix_t %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Convert back to a data frame 
  mean_expression_per_sample <- as.data.frame(mean_expression_per_sample)
  rownames(mean_expression_per_sample) <- mean_expression_per_sample$ID
  mean_expression_per_sample$ID <- NULL
  #head(t(mean_expression_per_sample))
  avg_expression_df <- t(mean_expression_per_sample)
  
  
  
  # Make design matrix. Need to make sure the order of the categories is the same
  # as the order of the columns in the average expression data frame
  # Get the correct order for the IR dose
  rejection_cat <- c()
  for (i in colnames(avg_expression_df)) {
    sample_id <- i # don't need to edit because sample ID was retained when pseudo bulking
    
    # get rejection type for this sample
    i_rt <- obj_CT@meta.data %>%
      filter(ID == sample_id) %>%
      dplyr::select(ID, IR.Dose) %>%
      unique() %>% dplyr::pull(IR.Dose)
    
    rejection_cat <- c(rejection_cat, i_rt)
  }
  rejection_cat <- as.factor(rejection_cat)
  
  if (as.numeric(table(rejection_cat)[1]) <3 || as.numeric(table(rejection_cat)[2]) <3) {
    print(paste("Skipping, ", ct, " because too few samples. Number of samples = ", length(rejection_cat), sep = ""))
  } else {
    # Set design matrix
    design <- model.matrix(~0+rejection_cat)
    
    # remove rejection_cat from column name to make cleaner
    colnames(design) <- gsub("rejection_cat", "", colnames(design))
    colnames(design) <- gsub(" ", "", colnames(design))
    colnames(design) <- gsub("20Gy", "Radiation", colnames(design))
    colnames(design) <- gsub("0Gy", "NoRadiation", colnames(design))
    
    # Set contrasts
    contr.matrix <- makeContrasts(
      noRadiationvsRadiation = Radiation-NoRadiation, # I think this will make upregulated in radiation have a pos logFC
      levels = colnames(design))
    #contr.matrix
    
    # Filtering low-exp genes  
    d0 <- DGEList(avg_expression_df)
    
    
    # I want to set a filter to keep genes that have expression in at least 3 samples
    keep <- rownames(avg_expression_df)[rowSums(avg_expression_df > 0) >= 3]
    #d <- d0[keep,] 
    keep_inter <- intersect(test_genes, keep)
    # remove hepatocyte genes from keep list
    # Genes to remove
    genes_to_remove <- c("Mup20", "Apoa2", "Serpina1e", "Arg1")
    
    filtered_genes <- setdiff(keep_inter, genes_to_remove)
    print(paste("Genes tested: ", length(filtered_genes), sep = ""))
    
    d <- d0[filtered_genes,] 
    
    dim(d0)
    dim(d)
    
    # NEW
    y <- voom(d, design, plot = T)
    fit <- lmFit(y, design)
    tmp <- contrasts.fit(fit, contr.matrix)
    tmp <- eBayes(tmp)
    
    #top.table <- topTable(tmp, sort.by = "P", n = Inf)
    
    # tmp is efit
    all_results[[ct]] <- topTable(tmp, number = 477)
    all_results[[ct]]$gene <- rownames(all_results[[ct]])
    all_results[[ct]]$cell_type <- ct
    all_results_efit[[ct]] <- tmp
  }
}


# Merge the list of data frames into one data frame
all_results_merged <- do.call(rbind, all_results)

all_results_merged_sig <- all_results_merged %>% filter(adj.P.Val <=0.05)

length(unique(all_results_merged_sig$cell_type))
length(unique(all_results_merged_sig$gene))
# 134 sig genes across 26 cell types


# Save results
#all_results <-readRDS("/scratch/aoill/projects/McGee_mouse_spatial/DGE_6_weeks_radiated_vs_not_new_filter.rds") 
saveRDS(all_results, "/scratch/aoill/projects/McGee_mouse_spatial/DGE_6_weeks_radiated_vs_not_new_filter.rds") 