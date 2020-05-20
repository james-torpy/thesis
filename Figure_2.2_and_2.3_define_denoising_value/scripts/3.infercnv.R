#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.2_define_denoising_value"
sample_name <- args[1]
numcores <- as.numeric(args[2])
analysis_mode <- args[3]
sim_name <- args[4]
denoise_value <- args[5]

#project_name <- "thesis"
#subproject_name <- "Figure_2.2_define_denoising_value"
#sample_name <- "CID4520N"
#numcores <- 6
#analysis_mode <- "samples"
#sim_name <- "sim7"
#denoise_value <- "1.3"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Number cores = ", numcores))
print(paste0("Analysis mode = ", analysis_mode))
print(paste0("Simulation name = ", sim_name))
print(paste0("Denoise value = ", denoise_value))

# if denoising, convert denoise_value into numeric class for infercnv:
if (denoise_value != "no") {
  denoise_value <- as.numeric(denoise_value)
}

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

out_path <- paste0(results_dir, "infercnv/", sample_name, "/", 
  sim_name, "/")
out_dir <- paste0(out_path, denoise_value, "_denoising/", 
  analysis_mode, "_mode/")
input_dir <- paste0(out_path, "/input_files/")

system(paste0("mkdir -p ", out_dir))

setwd(out_dir)

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))

print(paste0("Running InferCNV on ", sim_name))


################################################################################
### 1. Run InferCNV ###
################################################################################

# define normals which will act as InferCNV reference cells:
normals <- "Non_epithelial"
print(paste0("Normal is: ", normals))

print("Creating inferCNV object...")
raw_path <- paste0(input_dir, "input_matrix.txt")
annotation_path <- paste0(input_dir, "metadata.txt")
gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
initial_infercnv_object <- CreateInfercnvObject(
  raw_counts_matrix=raw_path,
  annotations_file=annotation_path,
  delim="\t",
  gene_order_file=gene_path,
  ref_group_names=normals
)
print("InferCNV object created, running inferCNV...")

if (denoise_value == "no") {
  system.time(
    infercnv_output <- try(
      infercnv::run(
        initial_infercnv_object,
        num_threads=numcores-1,
        out_dir=".",
        cutoff=0.1,
        window_length=101,
        max_centered_threshold=3,
        plot_steps=F,
        denoise=F,
        analysis_mode = analysis_mode
      )
    )
  )
  system("touch ./infercnv.png")
} else {
  system.time(
    infercnv_output <- try(
      infercnv::run(
        initial_infercnv_object,
        num_threads=numcores-1,
        out_dir=".",
        cutoff=0.1,
        window_length=101,
        max_centered_threshold=3,
        plot_steps=F,
        denoise=T,
        sd_amplifier=denoise_value,
        analysis_mode = analysis_mode
      )
    )
  )
}



