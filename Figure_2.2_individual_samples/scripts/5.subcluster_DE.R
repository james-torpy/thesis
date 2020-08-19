#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
specific_DE <- args[6]  # DE genes found for subclusters listed first
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
specific_features <- args[7]  
if (specific_features != "none") {
  specific_features <- strsplit(
    specific_features,
    "_"
  )[[1]]
}


#project_name <- "thesis"
#subproject_name <- "Figure_2.2_individual_samples"
#sample_name <- "CID3586"
#subcluster_method <- "random_trees"
#subcluster_p <- "0.05"
#if (subcluster_p != "none") {
#  subcluster_p <- as.numeric(subcluster_p)
#}
#coverage_filter <- "filtered"
#remove_artefacts <- "artefacts_not_removed"
#specific_DE <- "none"
##specific_DE <- "CNV_1.CNV_2..CNV_3.CNV_4.CNV_5.CNV_6"
#if (specific_DE != "none") {
#  specific_DE <- strsplit(
#    strsplit(
#      specific_DE,
#      "\\.\\."
#    )[[1]],
#    "\\."
#  )
#}
#specific_features <- "none"
##specific_features <- paste(
##  c(
##    "ELF5_GATA3_MDK_CXCL17_MGST1_TACSTD2_MUCL1_S100A14_S100A16_CTSD_ASCL1",
##    "SLC9A3R1_TFPI2_TMEM176A",
##    "IFI6_AREG_IFITM1_ZFP36"
##  ), collapse = "_"
##)
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
### 0. Define functions ###
################################################################################

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]


################################################################################
### 1. Load and format data ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(in_dir, "04_seurat_object_annotated.Rdata"))

# load metadata and subset:
epi_meta <- readRDS(
  paste0(
    Robject_dir, 
    "4b.final_epithelial_metadata_without_normals.Rdata"
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

if (!file.exists(paste0(table_dir, "subpop_DE.txt"))) {

  all_DE <- FindAllMarkers(
    only.pos = T,
    object = seurat_sub,
    min.pct = 0.5, 
    logfc.threshold = 0.75, 
    test.use = 'MAST'
  )

  all_DE_sorted <- arrange(all_DE, cluster, desc(avg_logFC))

  write.table(
    all_DE_sorted, 
    paste0(table_dir, "subpop_DE.txt"),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )

} else {

  all_DE_sorted <- read.table(
    paste0(table_dir, "subpop_DE.txt"),
    header = T
  )

}

heatmap_genes <- all_DE_sorted %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

if (!file.exists(paste0(plot_dir, "top_subpop_DE_heatmap.png"))) {
  hmap <- DoHeatmap(
    seurat_sub,
    features = heatmap_genes$gene,
    group.by = "ident"
  ) + theme(
    text = element_text(size = 20)
  )
  
  png(
    paste0(plot_dir, "top_subpop_DE_heatmap.png"),
    height = 9,
    width = 15,
    res = 300,
    units = "in"
  )
    print(hmap)
  dev.off()
}


################################################################################
### 2. DE between specific subclusters and plot ###
################################################################################

if (specific_DE != "none") {

  filename <- paste0(
    "DE_",
    paste(specific_DE[[1]], collapse = "_"),
    "_vs_",
    paste(specific_DE[[2]], collapse = "_")
  )
  
  if (!file.exists(paste0(plot_dir, filename))) {
  
    # merge Idents:
    custom_seurat <- seurat_sub
    custom_idents <- as.character(Idents(custom_seurat))
    custom_idents[
      custom_idents %in% specific_DE[[1]]
    ] <- "CNV population 1"
    custom_idents[
      custom_idents %in% specific_DE[[2]]
    ] <- "CNV population 2"
    Idents(custom_seurat) <- factor(custom_idents)
  
    # scale all genes:
    custom_seurat <- ScaleData(
      custom_seurat, 
      features = rownames(custom_seurat)
    )
  
    # find DE markers:
    custom_DE <- FindMarkers(
      object = custom_seurat,
      ident.1 = "CNV population 1",
      ident.2 = "CNV population 2",
      min.pct = 0.5, 
      logfc.threshold = 0.75, 
      test.use = 'MAST'
    )
    
    # split into up and downregulated genes and take top 20 of each:
    custom_DE_split <- split(
      custom_DE, 
      custom_DE$avg_logFC >0
    )
    custom_DEs_sorted <- lapply(custom_DE_split, function(x) {
      return(
        x %>%
          rownames_to_column(var = "gene") %>%
          arrange(desc(avg_logFC)) %>% 
          top_n(20, avg_logFC)
      )
    })
    custom_DE_sorted <- do.call("rbind", custom_DEs_sorted)

    write.table(
      custom_DE_sorted, 
      paste0(table_dir, "custom_pop_DE.txt"),
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  
    hmap <- DoHeatmap(
      custom_seurat,
      features = custom_DE_sorted$gene,
      group.by = "ident"
    ) + theme(
      text = element_text(size = 20)
    )
    
    png(
      paste0(plot_dir, "custom_pop_DE_heatmap.png"),
      height = 9,
      width = 15,
      res = 300,
      units = "in"
    )
      print(hmap)
    dev.off()
  
  }

  if (specific_features != "none") {

    # find DE markers:
    custom_gene_DE <- FindMarkers(
      object = custom_seurat,
      ident.1 = "CNV population 1",
      ident.2 = "CNV population 2",
      features = specific_features,
      logfc.threshold = 0.2, 
      test.use = 'MAST'
    )
  
    custom_gene_DE_sorted <- custom_gene_DE %>%
      rownames_to_column(var = "gene") %>%
      arrange(desc(avg_logFC)) %>%
      filter(p_val_adj < 0.1)
  
    write.table(
      custom_gene_DE_sorted, 
      paste0(table_dir, "custom_gene_pop_DE.txt"),
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  
    hmap <- DoHeatmap(
      custom_seurat,
      features = custom_gene_DE_sorted$gene,
      group.by = "ident"
    ) + theme(
      text = element_text(size = 20)
    )
    
    png(
      paste0(plot_dir, "custom_gene_DE_heatmap.png"),
      height = 9,
      width = 15,
      res = 300,
      units = "in"
    )
      print(hmap)
    dev.off()

  }

}



