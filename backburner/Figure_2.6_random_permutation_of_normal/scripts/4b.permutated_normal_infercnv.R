#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
sample_name <- args[1]
numcores <- as.numeric(args[2])
analysis_mode <- args[3]
permutation_proportion <- args[4]
simulation_number <- as.character(args[5])

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
sample_name <- "permutated_CID4520N"
numcores <- 10
analysis_mode <- "samples"
permutation_proportion <- 0.2
simulation_number <- "1"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Number cores = ", numcores))
print(paste0("Analysis mode = ", analysis_mode))
print(paste0("Permutation proportion = ", permutation_proportion))
print(paste0("Simulation number = ", simulation_number))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

out_path <- paste0(results_dir, "infercnv/",
  sample_name, "/", permutation_proportion, "_proportion/",
  simulation_number, "/"
)

out_dir <- paste0(out_path, "/", analysis_mode, "_mode/")
system(paste0("mkdir -p ", out_dir))
input_dir <- paste0(out_path, "input_files/")

setwd(out_dir)

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))

print(paste0("Running InferCNV on ", sample_name , " replicate ", 
  simulation_number, ", assessing normal breast dataset with ", 
  permutation_proportion, " of genes permutated..."))


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
if (exists("initial_infercnv_object")) {
  print("InferCNV object created, running inferCNV...")
} else {
  print("InferCNV object not created, check inputs...")
}

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

