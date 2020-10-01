#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))
numcores <- as.numeric(args[2])
print(paste0("Numcores = ", numcores))
analysis_mode <- args[3]
print(paste0("Analysis mode = ", analysis_mode))
denoise_value <- args[4]
print(paste0("Denoise value = ", denoise_value))

sample_name <- "CID4515"
numcores <- 6
nUMI_threshold <- 25000
nGene_threshold <- 5000
analysis_mode <- "samples"
denoise_value <- "no"

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

out_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  analysis_mode, "_mode/", denoise_value, "_denoising/")
Robject_dir <- paste0(out_dir, "Rdata/")
plot_dir <- paste0(out_dir, "plots/")

input_dir <- paste0(out_dir, "input_files/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", out_dir))
system(paste0("mkdir -p ", input_dir))

setwd(out_dir)

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name ,
  "..."))


################################################################################
### 0. Define functions ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))


################################################################################
### 1. Fetch counts matrix from normal sample and create original infercnv input 
# files  ###
################################################################################

if (
  !file.exists(paste0(input_dir, "/input_matrix.txt")) | 
  !file.exists(paste0(input_dir, "/metadata.txt"))
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

  ######
  count_df <- count_df[1:5000,]
  ######

  # create metadata df:
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, subset_data=subset_data, 
    count_df, for_infercnv=T)
  seurat_10X <- infercnv_metadata$seurat
  
  # only keep cells in metadata df:
  print("Dimensions of count df before filtering out CAFs/Unknown cells: ") 
  print(dim(count_df))
  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  print("Dimensions of count df after filtering out CAFs/Unknown cells: ") 
  print(dim(count_df))

  print("Dimensions of metadata df after filtering out CAFs/Unknown cells: ") 
  print(dim(infercnv_metadata$metadata))

  # save original metadata:
  saveRDS(infercnv_metadata, paste0(Robject_dir, "/original_infercnv_metadata.Rdata"))
  
  # change non-epithelial cell type to "Non_epithelial":
  final_metadata <- infercnv_metadata$metadat
  final_metadata$cell_type[grep("pithelial", final_metadata$cell_type)] <- "Epithelial"
  final_metadata$cell_type[
    grep("pithelial", final_metadata$cell_type, invert=T)
  ] <- "Non_epithelial"

  # write infercnv input files:
  write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  write.table(final_metadata, paste0(input_dir, "metadata.txt"), 
    sep = "\t", quote = F, col.names = F, row.names = F)

}


################################################################################
### 2. Run InferCNV ###
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

  if (denoise_value != "no") {
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
          denoise=F,
          analysis_mode = analysis_mode
        )
      )
    )
  }
  

} else {
  print("InferCNV object not created properly, check and rerun...")
}


################################################################################
### 3. Plot distribution of all and neutral CNV signal ###
################################################################################

infercnv_output <- unlist(
  read.table(
    paste0(out_dir, "infercnv.observations.txt", 
    as.is=T, 
    header = T
  )
)
  
infercnv_df <- data.frame(score = infercnv)
p <- ggplot(infercnv_df, aes(score)) + 
  geom_density() + 
  geom_vline(xintercept = 0.91, colour = "red") + 
  geom_vline(xintercept = 1.09, colour = "red")
pdf(paste0(out_dir, "all_cancer_CNV_score_density_plot.pdf"))
  p
dev.off()
png(paste0(out_dir, "all_cancer_CNV_score_density_plot.png"))
  p
dev.off()

# isolate neutral peak and surrounds:
rough_neutral_df <- data.frame(
  score = infercnv_df[infercnv_df$score > 0.95 & infercnv_df$score < 1.05,]
)

# plot ad mark 1 and 2 std devs from mean:
p <- ggplot(rough_neutral_df, aes(score)) + 
  geom_density() 
#  + 
#  geom_vline(xintercept = 0.955, colour = "red") + 
#  geom_vline(xintercept = 1.045, colour = "red")
pdf(paste0(out_dir, "all_cancer_CNV_neutral_score_density_plot.pdf"))
  p
dev.off()
png(paste0(out_dir, "all_cancer_CNV_neutral_score_density_plot.png"))
  p
dev.off()


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

