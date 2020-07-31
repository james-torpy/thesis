#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# nUMI
# nGene
# Expression clusters
# CNV subclusters

#### Generate following t-SNEs/UMAPS: ###
# Best matching resolution clusters annotated by cell type
# Broad cell marker feature plot
# Epithelial vs myoepithelial feature plot
# Luminal vs basal feature plot

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
all_cell_PC <- args[2]
all_cell_res <- args[3]
malignant_PC <- args[4]
malignant_res <- args[5]
broad_markers <- unlist(
  strsplit(
    args[6],
    "_"
  )
)
epithelial_markers <- unlist(
  strsplit(
    args[7],
    "_"
  )
)
nUMI_threshold <- as.numeric(args[8])
nGene_threshold <- as.numeric(args[9])
min_proportion_cells_for_subcluster <- as.numeric(args[10])
random_tree_p <- as.numeric(args[11])
subcluster_by <- args[12]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("All cell PC = ", all_cell_PC))
print(paste0("All cell res = ", all_cell_res))
print(paste0("Malignant PC = ", malignant_PC))
print(paste0("Malignant res = ", malignant_res))
print("Broad markers: ")
print(broad_markers)
print("Epithelial markers: ")
print(epithelial_markers)
print(paste0("nUMI threshold = ", nUMI_threshold))
print(paste0("nGene threshold = ", nGene_threshold))
print(paste0("Min proportion of cells for subcluster to be defined = ", 
  min_proportion_cells_for_subcluster))
print(paste0("Subclustering by ", subcluster_by))

sample_name <- "CID4515"
all_cell_PC <- "C"
all_cell_res <- "PC_C_res.0.6"
malignant_PC <- "D"
malignant_res <- "SUBSET_D_res.0.8"
broad_markers <- unlist(
  strsplit(
    c(
      "ACTB_PTPRC_CD19_CD3D_CD68_PDGFRB_PECAM1_EPCAM"
    ),
    "_"
  )
)
epithelial_markers <- unlist(
  strsplit(
    c(
      "ACTB_EPCAM_KRT18_KRT5_MKI67"
    ),
    "_"
  )
)
nUMI_threshold <- as.numeric("8000")
nGene_threshold <- as.numeric("1300")
min_proportion_cells_for_subcluster <- as.numeric("0.005")
random_tree_p <- 0.05
subcluster_by <- "random_trees"

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(MAST, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(ggplot2)
library(Seurat)
library(infercnv, lib.loc = lib_loc)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(fpc, lib.loc = lib_loc)
library(dplyr)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
in_dir <- paste0(results_dir, "infercnv/", sample_name, "/p_", 
  random_tree_p, "/")
input_dir <- paste0(in_dir, "/input_files/")

out_dir <- paste0(in_dir, "clustered_by_", subcluster_by, "/")
Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]


#################################################################################
### 1. Load all cells seurat object and create tSNEs/UMAPs ###
#################################################################################

if (!file.exists(paste0(plot_dir, "epithelial_feature.png"))) {

  seurat_10X <- readRDS(paste0(seurat_dir, "03_seurat_object_processed.Rdata"))
  Idents(seurat_10X) <- paste0(
    seurat_10X@meta.data$garnett_seurat_cluster_call_major_PC_C_res.0.6,
    " ",
    seurat_10X@meta.data$PC_C_res.0.6
  )
  
  tSNE <- DimPlot(
    object = seurat_10X,
    label.size = 70,
    pt.size = 2,
    reduction = paste0("TSNE", all_cell_PC),
     order = T
  )
  png(
    paste0(plot_dir, all_cell_res, "_tsne.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(tSNE)
  dev.off()
  
  UMAP <- DimPlot(
    object = seurat_10X,
    label.size = 70,
    pt.size = 2,
    reduction = paste0("UMAP", all_cell_PC),
    order = T
  )
  
  png(
    paste0(plot_dir, all_cell_res, "_umap.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(UMAP)
  dev.off()
  
  broad_feature <- FeaturePlot(
    object = seurat_10X,
    features = broad_markers,
    pt.size = 1.5,
    reduction = paste0("TSNE", all_cell_PC),
    order = T
  )
  dev.off()
  png(
    paste0(plot_dir, "broad_feature.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(broad_feature)
  dev.off()
  
  epithelial_feature <- FeaturePlot(
    object = seurat_10X,
    features = epithelial_markers,
    pt.size = 1.5,
    reduction = paste0("TSNE", all_cell_PC),
    order = T
  )
  dev.off()
  png(
    paste0(plot_dir, "epithelial_feature.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(epithelial_feature)
  dev.off()

}


#################################################################################
### 1. Load malignant epithelial seurat object and create tSNEs/UMAPs to choose 
# best PC/resolution ###
#################################################################################

if (file.exists(paste0(seurat_dir, "05_seurat_object_malignant.Rdata"))) {

  seurat_malignant <- readRDS(paste0(seurat_dir, "05_seurat_object_malignant.Rdata"))
  Idents(seurat_malignant) <- paste0(
    seurat_malignant@meta.data$garnett_call_ext_major,
    " ",
    eval(parse(text = paste0("seurat_malignant@meta.data$", malignant_res)))
  )
  levels(Idents(seurat_malignant)) <- naturalsort(levels(Idents(seurat_malignant)))
  
  # define cluster annotation colours:
  cluster_number <- length(unique(Idents(seurat_malignant)))
  cluster_cols <- col_palette[1:cluster_number]
  names(cluster_cols) <- unique(Idents(seurat_malignant))

  saveRDS(cluster_cols, paste0(Robject_dir, "cluster_colours.Rdata"))
  
  if (!file.exists(paste0(plot_dir, malignant_res, "_epithelial_tsne.png"))) {
  
    # isolate and plot epithelial cells:
    epi_tSNE <- DimPlot(
      object = seurat_malignant,
      label.size = 40,
      pt.size = 2,
      reduction = paste0("TSNE", malignant_PC),
      cols = cluster_cols
    )
    
    png(
      paste0(plot_dir, malignant_res, "_epithelial_tsne.png"),
       width = 10, 
       height = 8, 
       res = 300, 
       units = 'in'
    )
      print(epi_tSNE)
    dev.off()
    
    epi_UMAP <- DimPlot(
      object = seurat_malignant,
      label.size = 40,
      pt.size = 2,
      reduction = paste0("UMAP", malignant_PC),
      cols = cluster_cols
    )
    
    png(
      paste0(plot_dir, malignant_res, "_epithelial_umap.png"),
       width = 10, 
       height = 8, 
       res = 300, 
       units = 'in'
    )
      print(epi_UMAP)
    dev.off()
  
  }
  
#  if (!file.exists(paste0(plot_dir, "expression_cluster_DGE_heatmap.png"))) {
#  
#    epi_DE <- FindAllMarkers(
#      only.pos = T,
#      object = seurat_malignant,
#      min.pct = 0.5, 
#      logfc.threshold = 0.7, 
#      test.use = 'MAST'
#    )
#    
#    epi_DE_sorted <- arrange(epi_DE, cluster, desc(avg_logFC))
#    
#    write.table(epi_DE_sorted, paste0(table_dir, "expression_cluster_DGE.txt"))
#    
#    heatmap_genes <- epi_DE_sorted %>% 
#      group_by(cluster) %>% 
#      top_n(10, avg_logFC)
#    
#    png(
#      paste0(plot_dir, "expression_cluster_DGE_heatmap.png"),
#      width = 12,
#      height = 9,
#      res = 300,
#      units = "in"
#    )
#      print(DoHeatmap(
#        seurat_malignant,
#        features = heatmap_genes$gene,
#        group.by = "ident",
#        group.colors = cluster_cols,
#        size = 3
#      ))
#    dev.off()
#    
#  }

}