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
normals_removed <- as.logical(args[6])
epi_res <- args[7]
epi_PC <- args[8]
garnett_slot <- args[9]
remove_outliers <- as.logical(args[10])

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- "CID4463"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
normals_removed <- TRUE
epi_res <- "PC_C_res.1"
epi_PC <-"C"
remove_outliers <- FALSE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(naturalsort, lib.loc = lib_loc)
library(Seurat)
library(ggplot2)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
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
### 0. Define colours ###
################################################################################

expr_cols <- read.table(
  paste0(ref_dir, "expression_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

################################################################################
### 1. Load and format data ###
################################################################################

print("Loading metadata df...")

if (normals_removed) {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/4b.final_epithelial_metadata_without_normals.Rdata"))
} else {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/4b.final_epithelial_metadata_with_normals.Rdata"))
}
subcluster_meta <- subset(
  epithelial_metadata, 
  select = c(cell_ids, subcluster_id)
)

# load seurat object:
seurat_epi <- readRDS(
  paste0(
    seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata"
  )
)

# choose resolution if needed:
if (epi_res != "none") {
  Idents(seurat_epi) <- eval(
    parse(
      text = paste0("seurat_epi@meta.data$", epi_res)
    )
  )
}


################################################################################
### 2. Plot UMAP of reclustered epithelial cells ###
################################################################################

# order idents so epithelial cells are printed on top:
ident_order <- levels(Idents(seurat_epi))

if (remove_outliers) {
  # remove outliers - defined as those > 5 standard devs from mean:
  embeddings <- eval(parse(
    text = paste0("seurat_epi@reductions$UMAP", epi_PC, "@cell.embeddings")
  ))
  mean_embeddings <- apply(embeddings, 2, mean)
  std_dev_embeddings <- apply(embeddings, 2,sd)
  
  outliers1 <- rownames(embeddings)[
    embeddings[,1] > (mean_embeddings[1] + (5*std_dev_embeddings[1])) | 
    embeddings[,1] < (mean_embeddings[1] - (5*std_dev_embeddings[1]))
  ]
  if (length(outliers1) < 5) {
    include_cells <- rownames(embeddings)[
      !(rownames(embeddings) %in% outliers1)
    ]
  }
  
  outliers2 <- rownames(embeddings)[
    embeddings[,2] > (mean_embeddings[2] + (5*std_dev_embeddings[2])) | 
    embeddings[,2] < (mean_embeddings[2] - (5*std_dev_embeddings[2]))
  ]
  if (length(outliers2) < 5) {
    include_cells <- include_cells[
      !(include_cells %in% outliers2)
    ]
  }
  
  # remvoe epithelial clusters with <10 cells:
  split_idents <- split(
    Idents(seurat_epi),
    Idents(seurat_epi)
  )
  small_clusters <- levels(Idents(seurat_epi))[
    unlist(
      lapply(split_idents, function(x) length(x) < 10)
    )
  ]
  include_cells <- include_cells[
    include_cells %in% names(Idents(seurat_epi))[
      !(Idents(seurat_epi) %in% small_clusters)
    ]
  ]
} else {
  include_cells <- rownames(seurat_epi@meta.data)
}

include_cells <- rownames(seurat_epi@meta.data)[
  rownames(seurat_epi@meta.data) %in% subcluster_meta$cell_ids
]

epi_umap <- DimPlot(
  seurat_epi,
  cells = include_cells,
  pt.size = 1.5,
  reduction = paste0("UMAP", epi_PC),
  label = F
) + 
scale_color_manual(
  labels = paste0("Epithelial ", 1:length(ident_order)), 
  values = expr_cols[1:length(ident_order)]
)

png(
  file = paste0(plot_dir, "epithelial_cell_expr_clusters_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(epi_umap)
dev.off()


################################################################################
### 2. Plot UMAP owith CNV subclusters labelled ###
################################################################################

# update idents with subcluster ids:
seurat_epi@meta.data$subcluster_ids <- NA
seurat_epi@meta.data[subcluster_meta$cell_ids,]$subcluster_ids <- 
  as.character(subcluster_meta$subcluster_id)

Idents(seurat_epi) <- factor(
  as.character(seurat_epi@meta.data$subcluster_ids),
  levels = naturalsort(unique(as.character(seurat_epi@meta.data$subcluster_ids)))
)

labs <- gsub(
  "_", " ", unique(
    Idents(seurat_epi)[
      !is.na(Idents(seurat_epi))
    ]
  )
)

epi_umap <- DimPlot(
  seurat_epi,
  cells = include_cells,
  pt.size = 1.5,
  reduction = paste0("UMAP", epi_PC),
  label = F
) + 
scale_color_manual(
  labels = labs, 
  values = subcluster_cols[1:length(labs)]
)

png(
  file = paste0(plot_dir, "epithelial_cell_CNV_subclusters_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(epi_umap)
dev.off()




