#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.2_define_denoising_value"
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))
numcores <- as.numeric(args[2])
print(paste0("Numcores = ", numcores))
analysis_mode <- args[3]
print(paste0("Analysis mode = ", sample_mode))
denoise_value <- as.numeric(args[4])
print(paste0("Denoise value = ", denoise_value))

sample_name <- "PDX4386"
numcores <- 6
analysis_mode <- "samples"
denoise_value <- "1.5"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Number cores = ", numcores))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(ggplot2)
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
results_dir <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")

in_dir <- paste0(project_dir, "raw_files/seurat_objects/", sample_name, "/")

out_path <- paste0(results_dir, "infercnv/", sample_name, "/")
Robject_dir <- paste0(out_path, "/Rdata/")
plot_dir <- paste0(out_path, "/plots/")
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

input_dir <- paste0(out_path, "input_files/")
system(paste0("mkdir -p ", input_dir))

out_dir <- paste0(out_path, analysis_mode, "_mode/", denoise_value, "_denoising/")
system(paste0("mkdir -p ", out_dir))

setwd(out_dir)

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("InferCNV input file directory = ", input_dir))
print(paste0("InferCNV output file directory = ", out_dir))


print(paste0("Preparing normal breast InferCNV input files for ", sample_name,
  "..."))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))


################################################################################
### 1. Fetch counts matrix from normal sample and create infercnv input 
# files ###
################################################################################

if (
  !file.exists(paste0(input_dir, "input_matrix.txt")) | 
  !file.exists(paste0(input_dir, "metadata.txt"))
) {

  # load seurat object:
  seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
  
  # create raw matrix input file and subset if necessary:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
  print(
    paste0(
      "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
    )
  )

  count_df <- count_df[1:6000,]

  # create metadata df:
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, count_df, for_infercnv=T, 
    garnett="garnett_call_ext")
  seurat_10X <- infercnv_metadata$seurat

  cell_types <- paste(unique(infercnv_metadata$metadata$cell_type), collapse = " ")
  print(paste0("Cell types are: ", cell_types))
 
  # only keep cells in metadata df:
  print(paste0("No cells in count df before filtering for those in metadata df = ", 
      ncol(count_df)))
  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  print(paste0("No cells in count df after filtering for those in metadata df = ", 
      ncol(count_df)))

  # change non-epithelial cell type to "Non_epithelial":
  infercnv_metadata$metadata <- infercnv_metadata$metadat
  infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)] <- "Epithelial"
  infercnv_metadata$metadata$cell_type[
    grep("pithelial", infercnv_metadata$metadata$cell_type, invert=T)
  ] <- "Non_epithelial"

  # save normal dataset infercnv input files:
  write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
    sep = "\t", quote = F, col.names = F, row.names = F)

}


################################################################################
### 2. Run denoised InferCNV ###
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
        analysis_mode=analysis_mode
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
        sd_amplifier=denoise_value,
        analysis_mode=analysis_mode
      )
    )
  )
}


