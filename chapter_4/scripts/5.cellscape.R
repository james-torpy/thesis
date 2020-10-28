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
library(naturalsort, lib.loc = lib_loc)
library("cellscape", lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
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

add_chromosome_indices <- function(indices_df, chr_positions) {
  indices_df$chr_start <- chr_positions$index[indices_df$start]
  indices_df$chr_end <- chr_positions$index[indices_df$end]
  return(indices_df)
}


################################################################################
### 1. Load CNV and chromosome data and format ###
################################################################################

chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))

subpop_CNV_only <- readRDS(
  paste0(Robject_dir, "CNV_indices_and_lengths.Rdata")
)

epithelial_heatmap <- readRDS(
  paste0(Robject_dir, "5a.final_epithelial_heatmap_without_normals.Rdata")
)

# annotate all genes with chromosome position:
for (l in 1:length(chr_data$lengths)) {
  if (l==1) {
    chrs <- rep(names(chr_data$lengths)[l], chr_data$lengths[l])
  } else {
    chrs <- c(chrs, rep(names(chr_data$lengths)[l], chr_data$lengths[l]))
  }
}
chr_positions <- data.frame(
  gene = colnames(epithelial_heatmap),
  chr = chrs
)
chr_positions <- do.call(
  "rbind",
  lapply(
    split(chr_positions, chr_positions$chr),
    function(x) {
      x$index <- 1:nrow(x)
      return(x)
    }
  )
)
chr_positions <- chr_positions[naturalorder(chr_positions$chr),]
chr_positions$overall_index <- 1:nrow(chr_positions)

# add chromosomal indices:
subpop_CNV_only <- lapply(subpop_CNV_only, add_chromosome_indices, chr_positions)


  