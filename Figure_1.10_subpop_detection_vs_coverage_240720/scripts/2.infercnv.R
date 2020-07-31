#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_1.9_subpop_detection_vs_coverage"
sample_name <- args[1]
numcores <- as.numeric(args[2])
simulation_number <- args[3]
analysis_mode <- args[4]
p_val <- as.numeric(args[5])
downsample_proportion <- as.character(args[6])

project_name <- "thesis"
subproject_name <- "Figure_1.9_subpop_detection_vs_coverage"
sample_name <- "CID4520N_cancer_sim"
numcores <- 20
simulation_number <- "1"
analysis_mode <- "subclusters"
p_val <- 0.1
downsample_proportion <- "no"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Number cores = ", numcores))
print(paste0("Simulation number = ", simulation_number))
print(paste0("Analysis mode = ", analysis_mode))
print(paste0("p-value = ", p_val))
print(paste0("Downsample proportion = ", downsample_proportion))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/", simulation_number, "/", downsample_proportion, "_downsampling/")
input_dir <- paste0(in_dir, "/input_files/")

setwd(in_dir)

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("In/Output directory = ", in_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name ,
  " number ", simulation_number, ", downsampled to ", 
  downsample_proportion))


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

if (analysis_mode == "subclusters") {
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
        sd_amplifier=1.3,
        analysis_mode = analysis_mode,
        tumor_subcluster_pval = p_val
      )
    )
  )
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
        sd_amplifier=1.3,
        analysis_mode = analysis_mode
      )
    )
  )  
}

## remove temp files:
#system(paste0("rm ", out_dir, "/*dat"))
#system(paste0("rm ", out_dir, "/*preliminary*"))
#system(paste0("rm ", out_dir, "/01_*"))
#system(paste0("rm ", out_dir, "/02_*"))
#system(paste0("rm ", out_dir, "/03_*"))
#system(paste0("rm ", out_dir, "/04_*"))
#system(paste0("rm ", out_dir, "/05_*"))
#system(paste0("rm ", out_dir, "/06_*"))
#system(paste0("rm ", out_dir, "/07_*"))
#system(paste0("rm ", out_dir, "/08_*"))
#system(paste0("rm ", out_dir, "/09_*"))
#system(paste0("rm ", out_dir, "/10_*"))
#system(paste0("rm ", out_dir, "/11_*"))
#system(paste0("rm ", out_dir, "/12_*"))

