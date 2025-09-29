Proximity Analysis 1.0 - HOW TO RUN

Proximity analysis functions: spatial_proximity_functions.R
Proximity analysis code: spatial_proximity_v1.R
Proximity analysis bash for array run: run_spatial_proximity.sh

The code is set to run as an array job on each individual cell 
Input required is a data frame with the following:
1. spatial coordinates 
2. cell labels
3. sample ids 
4. unique cell ids as rownames 
#                    x_centroid y_centroid             cell_id               final_CT   sample
#VUILD106_aaaaiafn-1   2620.048   6594.044 VUILD106_aaaaiafn-1                    AT2 VUILD106
#VUILD110_aaabapik-1   3787.849   4875.265 VUILD110_aaabapik-1 Activated Fibrotic FBs VUILD110
#VUILD115_aaabdnbl-1   3353.458   3944.272 VUILD115_aaabdnbl-1                B cells VUILD115
#VUILD106_aaaampcm-1   4893.292   2680.374 VUILD106_aaaampcm-1                B cells VUILD106
NOTE: column names are rigid for this analysis, find and replace variable names in both R code if other column names are desired


# Step 0 - edit input and output directories in spatial_proximity_v1.R
----------------------------------------------------------------------------------------------------------------------------------
In spatial_proximity_v1.R change the following
LINE 17: add input directory
inDir <- '/data/'
LINE 18: add output directory
outDir <- '/output/'
LINE 22: change meta data file name and path if other than outDir
x <- read.csv(file.path(inDir, "proximity_analysis_example.csv"), row.names = 1)


# Step 1 - get a list of cell types
----------------------------------------------------------------------------------------------------------------------------------
inDir <- '/data/'
meta = read.csv(file.path(inDir, "proximity_analysis_example.csv"), row.names = 1)
dt <- table(meta$final_CT) %>% as.data.frame()
write.table(dt[,1], file.path(inDir, "example_celltype_list.txt"),
            append = F, quote = F,row.names = F,col.names = F, sep = "\t")
            
            
# Step 2 - run proximity analysis on all cells 
----------------------------------------------------------------------------------------------------------------------------------
cell_type_list=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/data/example_celltype_list.txt
container=rstudio-4.3.0-4-with_modules.sif  # change to your container
proximity_functions_rscript=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/Spatial_proximity_functions.R
proximity_rscript=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/spatial_proximity_v1.R
radius=30 # radius of 30 is ~3 cells in each direction

# Get the number of cell types from the cell type list
num_lines=$(wc -l < $cell_type_list)
if [ "$num_lines" -eq 0 ]; then
  echo "Error: celltype list is empty; SEE README TO GENERATE"
  exit 1
fi

# Submit array job
job_id=$(sbatch --parsable --array=1-${num_lines}%4 \
  -p compute -n 26 --mem 200GB --nodes=1 -t 0-36:00 \
  -o "slurm.loop.%A_%a.out" -e "slurm.loop.%A_%a.err" \
  run_spatial_proximity.sh "$cell_type_list" "$container" "$proximity_rscript" "$proximity_functions_rscript" "$radius")


# step 3 - assess output
----------------------------------------------------------------------------------------------------------------------------------
The code will generate a single csv with proximity results for each cell type. Below is an example output
# sample output
# A tibble: 6 × 9
#  distance cell.a                cell.b                degree celltypeA celltypeB angle   idx sid       
#     <dbl> <chr>                 <chr>                  <dbl> <chr>     <chr>     <dbl> <int> <chr>     
#1     9.30 s2_GOLIATH_mpjfkfag-1 s2_GOLIATH_mpjfkmok-1   6.23 AT1       AT2          30     1 s2_BPD03  
#2    14.8  s2_GOLIATH_akaejale-1 s2_GOLIATH_gafnpjoj-1   6.52 AT1       AT2          30     1 s2_BPD03  
#3     6.39 s2_GOLIATH_ijgbgfch-1 s2_GOLIATH_ijgbcdoa-1  37.6  AT1       Arterial     30     1 s2_4moCT27
#4     6.12 s2_GOLIATH_hglbinnm-1 s2_GOLIATH_hglamlnh-1  33.6  AT1       Arterial     30     1 s2_4moCT27
distance - the distance in pixel or microns (depends on input)
cell.a - anchor cell id
cell.b - proximal cell
degree - the degree of which the cell.b is located relative to cell.a
celltypeA - annotation of cell.a
celltypeB - annotation of cell.b
angle - a categorical summary of a range of degrees (e.g 30 for all cells within 30-60 degrees)
idx - nth neighbor assignment for each angle category
sid - sample ID


# Comments
----------------------------------------------------------------------------------------------------------------------------------
To run on transcripts cell id would be transcript id and annotations will be gene names
