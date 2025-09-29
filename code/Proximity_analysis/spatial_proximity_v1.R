## proximity code clean

#if(!require(sf)){
#  install.packages("sf", lib = "/home/hnatri/R/rstudio-4.3.0-4-with_modules.sif/libs/")
#  library(sf)
#}

library(sf, lib.loc = "/home/hnatri/R/rstudio-4.3.0-4-with_modules.sif/libs/")
library(dplyr)
library(data.table)
library(doParallel)
library(tidyr)
library(splancs, lib.loc = "/home/hnatri/R/rstudio-4.3.0-4-with_modules.sif/libs/")

# For debugging

#container=/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif  # change to your container
#cell_type_list=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/data/example_celltype_list.txt
#proximity_functions_rscript=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/Spatial_proximity_functions.R
#proximity_rscript=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/spatial_proximity_v1.R
#radius=30 # radius of 30 is ~3 cells in each direction
# run proximity analysis for each cell type
#cell_type=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$cell_type_list")
#singularity exec $container Rscript $proximity_rscript "$cell_type" "$proximity_functions_rscript" "$radius"
#args <- c("B1",
#          "/home/hnatri/PIPAC_spatial/code/Proximity_analysis/Spatial_proximity_functions.R",
#          30)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No cell type provided. Usage: Rscript spatial_proximity.R <argument>", call. = FALSE)
}

## load data and define code paramters
#source(args[2])
source(args[2])
inDir <- '/home/hnatri/PIPAC_spatial/code/Proximity_analysis/data/'
outDir <- '/home/hnatri/PIPAC_spatial/code/Proximity_analysis/output/'
r = as.numeric(args[3])  ## set radius 30 is ~3 cells in each direction

# read object metadata
x <- read.csv(file.path(inDir, "PIPAC_proximity_input.csv"), row.names = 1)
head(x)
#                       x_centroid y_centroid               cell_id    final_CT        sample
# s2_GOLIATH_gifcdhkj-1 10178.5625  10266.350 s2_GOLIATH_gifcdhkj-1         AT1   s2_BL092222
# s2_GOLIATH_cdppbhno-1   615.7618   8388.647 s2_GOLIATH_cdppbhno-1 Fibroblasts s2_BL091922.B
# s2_GOLIATH_dihjbmdf-1   539.1301   8757.209 s2_GOLIATH_dihjbmdf-1 Fibroblasts s2_BL091922.B

# get sample IDs
sample_ids <- x %>% 
  dplyr::pull(sample) %>% as.character %>% unique() 

# code is set to run on an individual cell type 
# set cell type
cell_type_id <- args[1]
# cell_type_id <- 'AT1'
cell_type_id_out <- gsub("\\/", "\\_", cell_type_id)
print(cell_type_id)

# set up enviorment to run in parallel
nCores <- 1
cl <- makeCluster(nCores)
registerDoParallel(cores=nCores)


# Export functions, variables and packages to workers
cell_id = 0
cell_x = 0
cell_y = 0
clusterExport(cl = cl, 
              varlist = c("circleFun", "calc_d", "getDir", "assign.angle", "anchorPoint", 
                          "get_neighbors", "cell_id", "cell_x", "cell_y", "r"),
              envir = environment())
clusterEvalQ(cl = cl, {
  library(sf, lib.loc = "/home/hnatri/R/rstudio-4.3.0-4-with_modules.sif/libs/")
  library(dplyr)
  library(data.table)
  library(tidyr)
})


## Loop over samples
proximal_cell_pop <- lapply(1:length(sample_ids), function(i){
  # do par needs packaged loaded for each worker - do not delete
  library(sf, lib.loc = "/home/hnatri/R/rstudio-4.3.0-4-with_modules.sif/libs/")
  library(dplyr)
  library(data.table)
  library(doParallel)
  library(tidyr)
  library(splancs, lib.loc = "/home/hnatri/R/rstudio-4.3.0-4-with_modules.sif/libs/")
  
  sid <- sample_ids[i]; message("finding neighbors for sample:", sid, "\n")
  obj.sid <- subset(x, sample == sid)
  
  # 1. get cell coordinates
  meta <- obj.sid
  cells_coord <- meta %>% 
    filter(final_CT == cell_type_id)
  
  ## skip is sample does not have the right cells
  if(nrow(cells_coord) < 5){
    return(NULL)
  }
  
  cell_counts <- c()
  # r <- 30 ## set radius ~30 is 3 cells 
  # print(r)
  for(j in 1:nrow(cells_coord)){
    # print(j) # used for diagnosis
    cell_id <- rownames(cells_coord)[j]
    split_cell.id <- strsplit(cell_id, ";") %>% unlist()
    cell_x <- cells_coord$x_centroid[j]
    cell_y <- cells_coord$y_centroid[j]
    
    tmp.d <- get_neighbors(cell_x_coords = cell_x,
                           cell_y_coords = cell_y,
                           r = r,
                           cell_coords_df = meta,
                           cell_id = cell_id)
    
    if(is.null(tmp.d)){next}
    tmp.d$cell.a <- cell_id
    tmp.d <- tmp.d[!duplicated(tmp.d),]
    
    # index rows by distance + angle
    tmp.d <- tmp.d %>%
      group_by(angle) %>%
      mutate(idx = rank(V1, ties.method = "min")) %>%
      ungroup()
    
    tmp.d$sid <- sid
    
    cell_counts <- rbind(cell_counts, tmp.d)
  }
  #cell_counts$sid <- sid
  
  return(cell_counts)
})

proximal_cell_pop_clean <- proximal_cell_pop[
  sapply(proximal_cell_pop, function(x) is.data.frame(x) && !is.null(x))
]
proximal_cell_pop_clean <- do.call("rbind", proximal_cell_pop_clean)
colnames(proximal_cell_pop_clean)[1] = 'distance'
stopImplicitCluster() 
out_fname <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
saveRDS(proximal_cell_pop_clean, file.path(outDir, out_fname)) 

message("\n","\n", "job run complete for cell type:", cell_type_id, "\n")
