#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### DGE between different CNV-driven subclusters ###

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]

sample_name <- "CID4515"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(scales, lib.loc = lib_loc)
library(ggplot2)
library(Seurat)
library(MAST, lib.loc = lib_loc)
library(dplyr)
library(RColorBrewer)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
#ref_dir <- paste0(project_dir, "/refs/")
#func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
in_dir <- paste0(results_dir, "infercnv/CID4515/Rdata/")

table_dir <- paste0(results_dir, "infercnv/CID4515/tables/")
plot_dir <- paste0(results_dir, "infercnv/CID4515/plots/")

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]

subcluster_cols <- rev(col_palette)[1:8]


################################################################################
### 1. Load seurat object and CNV subcluster data ###
################################################################################


epithelial_metadata <- readRDS(
  paste0(
  	in_dir, "2b.epithelial_metadata_with_CNV_subclusters.Rdata"
  )
)


################################################################################
### 2. Plot distributions of CNV lengths across all datasets ###
################################################################################




