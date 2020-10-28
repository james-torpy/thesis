#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
min_CNV_length <- as.numeric(args[6])
min_CNV_proportion <- as.numeric(args[7])
normals_removed <- as.logical(args[8])

project_name <- "thesis"
subproject_name <- "chapter_4"
sample_name <- "CID4463"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(ggplot2)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)
library(Seurat)
library(dplyr)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
raw_dir <- paste0(project_dir, "raw_files/seurat_objects/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)

Robject_dir <- paste0(in_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

subtype_cols <- read.table(
  paste0(ref_dir, "subtype_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)


################################################################################
### 1. Load and format data ###
################################################################################

print("Loading heatmap and metadata dfs...")

epithelial_metadata <- readRDS(paste0(Robject_dir, 
  "/5b.final_epithelial_metadata_without_normals.Rdata"))

print(paste0(
  "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
  identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
))

# load seurat object:
seurat_10X <- readRDS(
  paste0(raw_dir, sample_name, "/03_seurat_object_processed.Rdata")
)

# create df with CNA value and ki67 expression:
plot_df <- data.frame(
  cell = epithelial_metadata$cell_ids,
  CNA_value = epithelial_metadata$CNA_value
)
count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
ki67_df <- data.frame(
  cell = names(count_df[rownames(count_df) == "MKI67",]),
  count = count_df[rownames(count_df) == "MKI67",]
)
plot_df <- merge(plot_df, ki67_df, by="cell")


################################################################################
### 2. Plot CNA value vs ki67 expression ###
################################################################################

p <- ggplot(data = plot_df, aes(x=CNA_value, y=count))
p <- p + geom_point()

pdf(paste0(plot_dir, "CNA_value_vs_ki67.pdf"))
  print(p)
dev.off()





