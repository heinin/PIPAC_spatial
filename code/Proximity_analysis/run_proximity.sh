#!/bin/bash

#export SIMG_FILE_NAME=rstudio-4.3.0-4-with_modules.sif

container=/packages/containers/RStudio/rstudio-4.3.0-4-with_modules.sif  # change to your container
cell_type_list=/home/hnatri/PIPAC_spatial/code/Proximity_analysis/data/example_celltype_list.txt
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

echo "Submitted job array with ID: $job_id"