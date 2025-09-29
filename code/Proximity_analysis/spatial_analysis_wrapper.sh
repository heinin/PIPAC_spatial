#!/bin/bash

# define run variables - container, Rscript & cell type list path
container=/scratch/nhadad/containers/geospatial_latest.sif
cell_type_list=/scratch/nhadad/lung_datasets/prox_analysis_av/code_share/data/example_celltype_list.txt
proximity_functions_rscript=/scratch/nhadad/lung_datasets/prox_analysis_av/code_share/Spatial_proximity_functions.R
proximity_rscript=/scratch/nhadad/lung_datasets/prox_analysis_av/code_share/spatial_proximity_v1.R
radius=30 # radius of 30 is ~3 cells in each direction

# Get the number of cell types from celltype_list.txt
num_lines=$(wc -l < $cell_type_list)
if [ "$num_lines" -eq 0 ]; then
  echo "Error: celltype list is empty; SEE README TO GENERATE"
  exit 1
fi


# Submit the job with dynamic array
job_id=$(sbatch --parsable --array=1-${num_lines}%4 \
  -p compute -n 26 --mem 300 --nodes=1 -t 0-36:00 \
  -o "slurm.loop.%A_%a.out" -e "slurm.loop.%A_%a.err" \
  run_spatial_proximity.sh "$cell_type_list" "$container" "$proximity_rscript" "$proximity_functions_rscript" "$radius")

echo "Submitted job array with ID: $job_id"
