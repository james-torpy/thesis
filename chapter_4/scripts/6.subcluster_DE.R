#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

### This script performs DE between each CNV-based subpopulation and plots top DE results, and top CNA-associated
# results ###

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
adj_p_cutoff <- as.numeric(args[6])
specific_DE <- args[7]  # DE genes found for subclusters listed first
# compared to those listed second
if (specific_DE != "none") {
  strsplit(
    strsplit(
      specific_DE,
      "\\.\\."
    )[[1]],
    "\\."
  )
}
specific_features <- args[8]  
if (specific_features != "none") {
  specific_features <- strsplit(
    specific_features,
    "_"
  )[[1]]
}

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
adj_p_cutoff <- as.numeric("0.1")
#specific_DE <- "none"
specific_DE <- "CNV_3.CNV_4.CNV_5.CNV_6..CNV_1.CNV_2"
if (specific_DE != "none") {
  specific_DE <- strsplit(
    strsplit(
      specific_DE,
      "\\.\\."
    )[[1]],
    "\\."
  )
  names(specific_DE) <- c("ELF5", "non-ELF5")
}
specific_features <- "none"

#specific_features <- paste(
# c(
#    "MUCL1_IGFBP5_NDRG1_ELF5_ELF3_MDK_CXCL14_LY6D_CCND1_DUSP1_TIMP1_SERPINF1_SERPINB4_S100A6_S100A14_S100A16"
# ), collapse = "_"
#)
#if (specific_features != "none") {
#  specific_features <- strsplit(
#    specific_features,
#    "_"
#  )[[1]]
#}

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))
print(paste0("Specific grouping for DE? ", specific_DE))
print(paste0("Specific genes for DE? ", specific_features))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(RColorBrewer)
library(naturalsort, lib.loc = lib_loc)
library(Seurat)
library(MAST, lib.loc = lib_loc)
library(dplyr)
library(ggplot2)
library(searcher, lib.loc=lib_loc)
library(tibble)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/")
in_dir <- paste0(seurat_dir, sample_name, "/")
raw_dir <- paste0(project_dir, "raw_files/seurat_objects/")
results_dir <- paste0(project_dir, "results/")
out_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/DE/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In dir = ", in_dir))
print(paste0("Out dir = ", out_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define colours and functions ###
################################################################################

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

plot_DE <- dget(paste0(func_dir, "plot_DE.R"))


################################################################################
### 1. Load and format data ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(in_dir, "04_seurat_object_annotated.Rdata"))

# load metadata and subset:
epi_meta <- readRDS(
  paste0(
    Robject_dir, 
    "5b.final_epithelial_metadata_without_normals.Rdata"
  )
)

subcluster_meta <- subset(epi_meta, select = c("cell_ids", "subcluster_id"))

# change seurat object Idents to subcluster ids:
temp_idents <- as.character(Idents(seurat_10X))
names(temp_idents) <- as.character(names(Idents(seurat_10X)))
m <- match(subcluster_meta$cell_ids, names(temp_idents))
temp_idents[m] <- as.character(subcluster_meta$subcluster_id)
Idents(seurat_10X) <- factor(
  as.character(temp_idents),
  levels = naturalsort(unique(temp_idents))
)

# subset object keeping only subcluster cells:
seurat_sub <- subset(
  seurat_10X,
  idents = as.character(unique(subcluster_meta$subcluster_id))
)


################################################################################
### 2. DE between all subclusters and plot ###
################################################################################

# stringent DE with only top significant DE genes:
plot_DE(
  seurat_object = seurat_sub,
  min_pct = 0.5,
  logfc_thresh = 0.7,
  return_thresh = 0.01,
  only.pos = FALSE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_plot = TRUE,
  raw_dir,
  table_dir,
  Robject_dir,
  ref_dir,
  func_dir,
  file_prefix = "top_subpop_DE"
)

# loose DE with all reported CNA-assoc genes (to be used for common gene DE
# heatmap):
plot_DE(
  seurat_object = seurat_sub,
  min_pct = 0.1,
  logfc_thresh = 0,
  return_thresh = 1,
  only.pos = FALSE,
  plot_dir,
  filter_sig = FALSE,
  CNA_assoc_plot = TRUE,
  raw_dir,
  table_dir,
  Robject_dir,
  ref_dir,
  func_dir,
  file_prefix = "all_subpop_DE"
)

# stringent DE between specified groups:
plot_DE(
  seurat_object = seurat_sub,
  min_pct = 0.5,
  logfc_thresh = 0.7,
  return_thresh = 0.01,
  only.pos = FALSE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_plot = TRUE,
  raw_dir,
  table_dir,
  Robject_dir,
  ref_dir,
  func_dir,
  file_prefix = paste0(names(specific_DE)[2], "_vs_", names(specific_DE)[2]),
  specific_DE = specific_DE,
  group1 = names(specific_DE)[1],
  group2 = names(specific_DE)[2]
)


## loose DE for use with iDEA:
#plot_DE(
#  seurat_object = seurat_sub,
#  min_pct = 0,
#  logfc_thresh = 0,
#  return_thresh = 1,
#  only.pos = FALSE,
#  plot_dir,
#  filter_sig = FALSE,
#  CNA_assoc_only = FALSE,
#  raw_dir,
#  table_dir,
#  Robject_dir,
#  ref_dir,
#  file_prefix = "all_subpop_DE"
#)

