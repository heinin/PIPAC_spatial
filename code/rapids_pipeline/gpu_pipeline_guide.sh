# /tgen_labs/banovich/personal_space/hnatri/rapids_pipeline/rapids_cuda11_EDM_2024-1-9.sif

gpusrunit() {
NODES=$1
CPUS=$2
RAM=$3
srun -p gpu-a100 --gres=gpu:1 -n $NODES --cpus-per-task=${CPUS} --mem=${RAM}G --pty bash
}

gpusrunit 1 1 128 # Get a GPU node

module load singularity

singularity shell --nv /tgen_labs/banovich/personal_space/hnatri/rapids_pipeline/rapids_cuda11_EDM_2024-1-9.sif  # Drop into the container interactively

source activate rpy2

# Make sure you are in the directory containing rapids_scanpy_funcs.py and utils.py before you start jupyter

jupyter lab --ip 0.0.0.0 --port 6120 --no-browser # Don't forget to add .rc.tgen.org before the colon

# http://g-h-1-8-13.rc.tgen.org:6120/lab?token=44e2ac756510406f5f1cc2248c9b2543987bf94ae378d420
# http://g-h-1-9-31.rc.tgen.org:6120:lab?token=c5bb0149f21d1df93f4d04392d67d347ff08c244e2f579af

# Start in the notebook called seurat_to_anndata.ipynb. This should be using the rpy2 kernel.

# Next move to the clustering.ipynb notebook. This should be using the rapids kernel.

# After clustering we leave jupyter and move to R.

# Currently using anndata_to_seurat.R

#source activate sceasy # Activate the sceasy conda env
#
#R

# Run the following in commandline R. (Change file names accordingly)
# Note: loses some metadata!

#library(sceasy)
#library(Seurat)
#use_condaenv("sceasy")
#
#seurat_object = sceasy::convertFormat("/scratch/hnatri/PIPAC/merged_spatial_filtered_clustered_NC50_NN20_PC20_AnnData.h5ad", from="anndata", to="seurat", main_layer="data")
#saveRDS(seurat_object, "/scratch/hnatri/PIPAC/merged_spatial_filtered_clustered_NC50_NN20_PC20_Seurat.rds")