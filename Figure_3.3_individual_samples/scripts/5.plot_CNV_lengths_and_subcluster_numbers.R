#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
args = commandArgs(trailingOnly=TRUE)

sample_names <- c("CID4066", "CID4515", "CID3921", "CID4067", 
   "CID45171", "CID3941", "CID4523", "CID3948", 
   "CID4463", "CID4530N", "CID3963", "CID44041", 
   "CID4465", "CID4535")

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
#library(scales, lib.loc = lib_loc)
library(ggplot2)
#library(Seurat)
#library(dplyr)
library(RColorBrewer)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
#ref_dir <- paste0(project_dir, "/refs/")
#func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
#raw_dir <- paste0(project_dir, "raw_files/")
#seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
in_path <- paste0(results_dir, "infercnv/")
out_path <- paste0(in_path, "/CNV_length_and_subcluster_distribution/")

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))

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

for (i in 1:length(sample_names)) {

  epithelial_metadata <- readRDS(
    paste0(
    	in_path, sample_names[i], 
    	"/Rdata/4.subcluster_revised_epithelial_metadata.Rdata"
    )
  )

  CNV_data <- do.call(
  	"rbind",
  	readRDS(
      paste0(
        in_path, sample_names[i], 
       	"/Rdata/CNV_indices_and_lengths.Rdata"
      )
    )
  )
  CNV_data$sample <- sample_names[i]

  if (i==1) {

  	cluster_no <- data.frame(
  	  sample = sample_names[i],
  	  number = length(unique(epithelial_metadata$subcluster_id))
  	)
  	
  	CNV_lengths <- subset(
  	  CNV_data, select = c(sample, call, length, genomic_length)
  	)

  } else {

  	cluster_no <- rbind(
  	  cluster_no,
  	  data.frame(
    	sample = sample_names[i],
    	number = length(unique(epithelial_metadata$subcluster_id))
      )
    )

    CNV_lengths <- rbind(
      CNV_lengths,
      subset(
  	    CNV_data, select = c(sample, call, length, genomic_length)
  	  )
  	)

  }

}


################################################################################
### 2. Plot distributions of CNV lengths across all datasets ###
################################################################################

######
CNV_lengths <- CNV_lengths[CNV_lengths$length < 500,]
######

dplot <- density(CNV_lengths$length)
png(
  paste0(plot_dir, "CNV_length_distribution.png"),
  width = 7, height = 5, units = "in", res = 300
)
  plot(
    dplot,
    main = "",
    xlab = "CNV length (number of genes)"
  )
dev.off()







