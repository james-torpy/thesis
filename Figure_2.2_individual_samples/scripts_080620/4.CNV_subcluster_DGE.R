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

seurat_malignant <- readRDS(
  paste0(
  	seurat_dir, "05_seurat_object_malignant.Rdata"
  )
)

epithelial_metadata <- readRDS(
  paste0(
  	in_dir, "2b.epithelial_metadata_with_CNV_subclusters.Rdata"
  )
)

# order metadata as in seurat object:
epithelial_metadata <- epithelial_metadata[
  rownames(seurat_malignant@meta.data),
]
print(paste0(
  "Are epithelial_metadata rownames in the same order as seurat object?? ",
  identical(rownames(epithelial_metadata), rownames(seurat_malignant@meta.data))
))

# add CNV subcluster ids to seurat metadata:
seurat_malignant@meta.data$CNV_subcluster <- epithelial_metadata$subcluster_id


################################################################################
### 2. Perform DGE between all CNV subclusters ###
################################################################################

Idents(seurat_malignant) <- seurat_malignant@meta.data$CNV_subcluster

all_DE <- FindAllMarkers(
  only.pos = T,
  object = seurat_malignant,
  min.pct = 0.5, 
  logfc.threshold = 0.7, 
  test.use = 'MAST'
)

all_DE_sorted <- arrange(all_DE, cluster, desc(avg_logFC))

write.table(all_DE_sorted, paste0(table_dir, "CNV_subpop_DGE.txt"))

heatmap_genes <- all_DE_sorted %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

png(
  paste0(plot_dir, "CNV_subpop_DGE_heatmap.png"),
  width = 7,
  height = 5,
  res = 300,
  units = "in"
)
  print(DoHeatmap(
    seurat_malignant,
    features = heatmap_genes$gene,
    group.by = "ident",
    group.colors = subcluster_cols,
    size = 3
  ))
dev.off()

pdf(
  paste0(plot_dir, "CNV_subpop_DGE_heatmap.pdf"),
  width = 7,
  height = 5
)
  print(DoHeatmap(
    seurat_malignant,
    features = heatmap_genes$gene,
    group.by = "ident",
    group.colors = subcluster_cols,
    size = 3
  ))
dev.off()


