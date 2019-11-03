#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.1a_infercnv_method_schematic"
sample_name <- "CID4515"
numcores <- 6

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
raw_dir <- seurat_path <- paste0(project_dir, "raw_files/input_files/")
in_dir <- paste0(raw_dir, sample_name, "/") 
results_dir <- paste0(project_dir, "results/")
out_dir <- paste0(results_dir, "infercnv/", sample_name, "/")
system(paste0("mkdir -p ", out_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name))


################################################################################
### 1. Generate input matrix and metadata files ###
################################################################################

normals <- "Stromal"

print("Creating inferCNV object...")
raw_path <- paste0(in_dir, "input_matrix.txt")
annotation_path <- paste0(in_dir, "metadata.txt")
gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
initial_infercnv_object <- CreateInfercnvObject(
  raw_counts_matrix=raw_path,
  annotations_file=annotation_path,
  delim="\t",
  gene_order_file=gene_path,
  ref_group_names=normals
)

print("InferCNV object created, running inferCNV...")
system.time(
  infercnv::run(
    initial_infercnv_object,
    num_threads=numcores-1,
    out_dir=out_dir,
    cutoff=0.1,
    window_length=101,
    max_centered_threshold=3,
    cluster_by_groups=F,
    plot_steps=T,
    denoise=T,
    sd_amplifier=1.3,
    analysis_mode = "samples"
  )
)
