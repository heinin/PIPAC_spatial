#!/bin/bash
#SBATCH -p compute
##SBATCH -n 26
#SBATCH --mem 300
#SBATCH --nodes=1
#SBATCH -t 0-36:00
#SBATCH -o slurm.loop.%A_%a.out
#SBATCH -e slurm.loop.%A_%a.err

# Load required module
module load singularity

cell_type_list="$1"
container="$2"
proximity_rscript="$3"
proximity_functions_rscript="$4"
radius="$5"

cat "$cell_type_list"

# run proximity analysis for each cell type
cell_type=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$cell_type_list")
singularity exec $container Rscript $proximity_rscript "$cell_type" "$proximity_functions_rscript" "$radius"
